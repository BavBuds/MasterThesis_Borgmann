#!/usr/bin/env python3
"""
Assembly Polishing Pipeline
Sequentially applies:
 1. Racon (4 rounds) with original long reads
 2. Medaka (for ONT data)
 3. NextPolish with short Illumina reads (R1/R2)

Designed to work with assemblies from FLYE/CANU pipeline.
Includes robust resume functionality to continue from interrupted runs.

Author: Based on pipeline format by Max Borgmann
Date: 03.04.2025

Steps:
 1. For each input assembly:
    A) Run 4 rounds of Racon polishing with original long reads
    B) Run Medaka polishing (ONT only)
    C) Run NextPolish with short Illumina reads (paired-end)
 2. After each polishing step, optionally run BUSCO:
    - with selected lineage (e.g., mollusca_odb10)
    - Parse the JSON to track improvements
 3. Also run seqkit stats to get GC% and other basic stats
 4. Summarize improvements at each stage
"""

import os
import sys
import argparse
import subprocess
import logging
import shutil
import glob
import json
import pandas as pd
import hashlib
import tempfile
from datetime import datetime

# --------------------------------------------------------------------
# 1. Logging Setup
# --------------------------------------------------------------------
logger = logging.getLogger("assembly_polishing")
logger.setLevel(logging.DEBUG)

def setup_logging(log_dir, log_name="assembly_polishing.log"):
    os.makedirs(log_dir, exist_ok=True)
    log_fp = os.path.join(log_dir, log_name)

    fh = logging.FileHandler(log_fp)
    fh.setLevel(logging.DEBUG)
    fmt_file = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fh.setFormatter(fmt_file)

    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    fmt_console = logging.Formatter('%(message)s')
    ch.setFormatter(fmt_console)

    logger.addHandler(fh)
    logger.addHandler(ch)

    return log_fp

# --------------------------------------------------------------------
# 2. State Management for Resume Functionality
# --------------------------------------------------------------------
def state_file_path(output_dir):
    """Returns the path to the state file."""
    return os.path.join(output_dir, "polish_pipeline_state.json")

def get_file_hash(file_path):
    """Get SHA256 hash of a file for integrity verification."""
    if not os.path.exists(file_path):
        return None
    
    sha256_hash = hashlib.sha256()
    with open(file_path, "rb") as f:
        for byte_block in iter(lambda: f.read(4096), b""):
            sha256_hash.update(byte_block)
    return sha256_hash.hexdigest()

def read_state(output_dir):
    """
    Read the current pipeline state.
    Returns a dict with the structure:
    {
        "assembly_name": {
            "last_update": "timestamp",
            "tech": "pacbio-hifi|nano-raw|...",
            "original_file": "path/to/original.fasta",
            "polishing_steps": {
                "racon_1": {"status": "completed", "path": "...", "hash": "..."},
                "racon_2": {"status": "completed", "path": "...", "hash": "..."},
                "racon_3": {"status": "completed", "path": "...", "hash": "..."},
                "racon_4": {"status": "completed", "path": "...", "hash": "..."},
                "medaka": {"status": "completed", "path": "...", "hash": "..."},
                "nextpolish": {"status": "completed", "path": "...", "hash": "..."}
            },
            "evaluations": {
                "original_busco": {"status": "completed", "result": {...}},
                "racon_4_busco": {"status": "completed", "result": {...}},
                "medaka_busco": {"status": "completed", "result": {...}},
                "nextpolish_busco": {"status": "completed", "result": {...}},
                "original_stats": {"status": "completed", "result": {...}},
                "racon_4_stats": {"status": "completed", "result": {...}},
                "medaka_stats": {"status": "completed", "result": {...}},
                "nextpolish_stats": {"status": "completed", "result": {...}}
            }
        }
    }
    """
    state_path = state_file_path(output_dir)
    if not os.path.exists(state_path):
        return {}
    
    try:
        with open(state_path, 'r') as f:
            return json.load(f)
    except (json.JSONDecodeError, FileNotFoundError):
        logger.warning(f"State file corrupted or not found. Starting fresh state.")
        return {}

def write_state(state, output_dir):
    """Write the pipeline state to disk with atomic operations."""
    state_path = state_file_path(output_dir)
    os.makedirs(output_dir, exist_ok=True)
    
    # First write to a temporary file, then rename to ensure atomic operation
    with tempfile.NamedTemporaryFile('w', dir=output_dir, delete=False) as tmp:
        json.dump(state, tmp, indent=2)
        tmp_path = tmp.name
    
    # Atomic rename operation
    shutil.move(tmp_path, state_path)
    logger.debug(f"Updated pipeline state file: {state_path}")

def update_state(output_dir, assembly_name, step, status, **kwargs):
    """
    Update the pipeline state for a specific step.
    
    Args:
        output_dir: Output directory containing the state file
        assembly_name: Assembly name
        step: Step name (e.g., racon_1, medaka, nextpolish)
        status: Status of the step (completed, failed)
        **kwargs: Additional data to store (path, hash, result, etc.)
    """
    state = read_state(output_dir)
    
    # Ensure nested structure exists
    if assembly_name not in state:
        state[assembly_name] = {
            "last_update": datetime.now().isoformat(),
            "polishing_steps": {},
            "evaluations": {}
        }
    
    # Determine whether this is a polishing step or evaluation
    if step in ['racon_1', 'racon_2', 'racon_3', 'racon_4', 'medaka', 'nextpolish']:
        if "polishing_steps" not in state[assembly_name]:
            state[assembly_name]["polishing_steps"] = {}
        state[assembly_name]["polishing_steps"][step] = {"status": status, **kwargs}
    else:
        if "evaluations" not in state[assembly_name]:
            state[assembly_name]["evaluations"] = {}
        state[assembly_name]["evaluations"][step] = {"status": status, **kwargs}
    
    state[assembly_name]["last_update"] = datetime.now().isoformat()
    
    # Write the updated state
    write_state(state, output_dir)

def is_step_completed(output_dir, assembly_name, step, verify_file=None):
    """
    Check if a step has been completed successfully.
    If verify_file is provided, also check if the file exists and hash matches.
    """
    state = read_state(output_dir)
    
    # Check if assembly exists in state
    if assembly_name not in state:
        return False
    
    # Check if step is in polishing_steps or evaluations
    step_data = None
    if step in ['racon_1', 'racon_2', 'racon_3', 'racon_4', 'medaka', 'nextpolish']:
        if "polishing_steps" in state[assembly_name]:
            step_data = state[assembly_name]["polishing_steps"].get(step)
    else:
        if "evaluations" in state[assembly_name]:
            step_data = state[assembly_name]["evaluations"].get(step)
    
    # Check if step exists and was completed
    if not step_data or step_data.get("status") != "completed":
        return False
    
    # If a file path is given for verification
    if verify_file and "path" in step_data:
        stored_path = step_data["path"]
        
        # Check if the file exists
        if not os.path.exists(stored_path):
            logger.warning(f"File marked as completed but doesn't exist: {stored_path}")
            return False
        
        # If hash is stored, verify the hash
        if "hash" in step_data:
            current_hash = get_file_hash(stored_path)
            if current_hash != step_data["hash"]:
                logger.warning(f"File hash mismatch for {stored_path}")
                return False
    
    return True

def get_completed_step_data(output_dir, assembly_name, step):
    """Retrieve data for a completed step from the state file."""
    state = read_state(output_dir)
    if assembly_name not in state:
        return None
    
    if step in ['racon_1', 'racon_2', 'racon_3', 'racon_4', 'medaka', 'nextpolish']:
        if "polishing_steps" in state[assembly_name]:
            return state[assembly_name]["polishing_steps"].get(step)
    else:
        if "evaluations" in state[assembly_name]:
            return state[assembly_name]["evaluations"].get(step)
    
    return None

def update_assembly_info(output_dir, assembly_name, tech=None, original_file=None):
    """Update general information about the assembly."""
    state = read_state(output_dir)
    
    if assembly_name not in state:
        state[assembly_name] = {
            "last_update": datetime.now().isoformat(),
            "polishing_steps": {},
            "evaluations": {}
        }
    
    if tech:
        state[assembly_name]["tech"] = tech
    
    if original_file:
        state[assembly_name]["original_file"] = original_file
    
    state[assembly_name]["last_update"] = datetime.now().isoformat()
    
    write_state(state, output_dir)

def get_assembly_info(output_dir, assembly_name):
    """Get general information about the assembly."""
    state = read_state(output_dir)
    if assembly_name in state:
        return state[assembly_name]
    return None

# --------------------------------------------------------------------
# 3. Utility: run commands, parse JSON, get seqkit stats
# --------------------------------------------------------------------
def run_cmd(cmd, cwd=None, desc=None):
    """Run a shell command with logging."""
    if desc:
        logger.info(f"🔧 {desc}")
    logger.info(f"Running: {cmd}")
    try:
        result = subprocess.run(cmd, shell=True, check=True, cwd=cwd,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if result.stdout.strip():
            logger.debug(f"STDOUT: {result.stdout.strip()}")
        if result.stderr.strip():
            logger.debug(f"STDERR: {result.stderr.strip()}")
        return True, result.stdout
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed: {cmd}")
        logger.error(f"Return code: {e.returncode}")
        logger.error(f"STDOUT: {e.stdout}")
        logger.error(f"STDERR: {e.stderr}")
        return False, e.stderr

def parse_busco_json(json_file):
    """
    Parse BUSCO v5+ JSON to extract:
      - 'Complete percentage' -> complete_percent
      - 'Missing percentage'  -> missing_percent
      - 'Contigs N50'         -> n50
      - 'Number of contigs'   -> num_contigs
      - 'Total length'        -> total_length
    """
    if not os.path.exists(json_file):
        logger.warning(f"BUSCO JSON not found: {json_file}")
        return {}

    with open(json_file, 'r') as f:
        try:
            js = json.load(f)
        except Exception as e:
            logger.error(f"Error loading BUSCO JSON {json_file}: {e}")
            return {}

    results = js.get("results", {})
    metrics = js.get("metrics", {})

    def to_int(x):
        try:
            return int(x)
        except:
            return 0
    def to_float(x):
        try:
            return float(x)
        except:
            return 0.0

    data = {}
    data["complete_percent"] = to_float(results.get("Complete percentage", 0.0))
    data["missing_percent"]  = to_float(results.get("Missing percentage", 0.0))
    data["n50"]             = to_int(metrics.get("Contigs N50", 0))
    data["num_contigs"]     = to_int(metrics.get("Number of contigs", 0))
    data["total_length"]    = to_int(metrics.get("Total length", 0))

    return data

def run_seqkit_stats(fasta_file, out_file=None):
    """
    Run seqkit stats on a FASTA file to get basic sequence statistics.
    Returns a dictionary with statistics.
    """
    if not os.path.exists(fasta_file):
        logger.error(f"FASTA file not found: {fasta_file}")
        return {}
    
    if not out_file:
        out_file = f"{fasta_file}.seqkit_stats.txt"
    
    cmd = f"seqkit stats -a {fasta_file} -o {out_file}"
    ok, _ = run_cmd(cmd, desc="seqkit stats")
    if not ok or not os.path.exists(out_file):
        logger.error(f"Failed to run seqkit stats on {fasta_file}")
        return {}
    
    try:
        df = pd.read_csv(out_file, sep='\t', comment='#')
        if df.empty:
            logger.warning(f"seqkit stats produced empty output for {fasta_file}")
            return {}
        
        # Convert the first row to a dictionary
        stats = df.iloc[0].to_dict()
        
        # Convert values to appropriate types
        for key in ['num_seqs', 'min_len', 'avg_len', 'max_len', 'sum_len']:
            if key in stats:
                stats[key] = int(stats[key])
        
        for key in ['gc_content']:
            if key in stats:
                stats[key] = float(stats[key])
        
        return stats
    
    except Exception as e:
        logger.error(f"Error parsing seqkit stats for {fasta_file}: {e}")
        return {}

def run_busco(fasta_file, out_dir, lineage, threads, assembly_name, step, state_dir=None, force=False):
    """
    Runs BUSCO on the given FASTA file.
    
    Args:
        fasta_file: Path to the assembly FASTA file
        out_dir: Directory to store BUSCO results
        lineage: BUSCO lineage to use (e.g., mollusca_odb10)
        threads: Number of CPU threads
        assembly_name: Name of the assembly
        step: Current polishing step (e.g., original, racon_4, medaka, nextpolish)
        state_dir: Directory for state management
        force: Force re-run even if already completed
    
    Returns:
        Path to BUSCO JSON result file or None if failed
    """
    step_key = f"{step}_busco"
    
    # Check if BUSCO has already been run
    if not force and state_dir and is_step_completed(state_dir, assembly_name, step_key):
        step_data = get_completed_step_data(state_dir, assembly_name, step_key)
        if step_data and "path" in step_data and "result" in step_data:
            json_file = step_data["path"]
            if os.path.exists(json_file):
                logger.info(f"✅ Resuming: BUSCO already completed for {assembly_name}/{step}")
                return json_file
    
    os.makedirs(out_dir, exist_ok=True)
    out_name = f"{assembly_name}_{step}"
    
    # Check if the BUSCO output directory already exists and has completed JSON
    busco_subdir = os.path.join(out_dir, out_name)
    json_candidates = glob.glob(os.path.join(busco_subdir, "short_summary.*.json"))
    
    if json_candidates and not force:
        logger.info(f"✅ Resuming: Found existing BUSCO results for {assembly_name}/{step}")
        json_file = json_candidates[0]
    else:
        # Run BUSCO
        cmd = (
            f"busco -i {fasta_file} -o {out_name} -l {lineage} -m genome "
            f"--cpu {threads} --out_path {out_dir} --force"
        )
        ok, _ = run_cmd(cmd, desc=f"BUSCO for {assembly_name}/{step}")
        if not ok:
            if state_dir:
                update_state(state_dir, assembly_name, step_key, "failed", 
                             error="busco_failed")
            return None

        # Check for JSON after running
        json_candidates = glob.glob(os.path.join(busco_subdir, "short_summary.*.json"))
        if not json_candidates:
            logger.error(f"No BUSCO JSON found for {assembly_name}/{step}")
            if state_dir:
                update_state(state_dir, assembly_name, step_key, "failed", 
                             error="busco_missing_json")
            return None
        json_file = json_candidates[0]
    
    # Parse the BUSCO results
    busco_results = parse_busco_json(json_file)
    
    # Update state
    if state_dir:
        update_state(
            state_dir, assembly_name, step_key, "completed", 
            path=json_file, 
            result=busco_results,
            timestamp=datetime.now().isoformat()
        )
    
    return json_file

# --------------------------------------------------------------------
# 4. Racon Polishing
# --------------------------------------------------------------------
def run_racon_polishing(assembly_fasta, reads_file, output_dir, assembly_name, round_num, threads=16, state_dir=None, force=False):
    """
    Run a single round of Racon polishing.
    
    Args:
        assembly_fasta: Input assembly FASTA
        reads_file: Long reads file (FASTQ)
        output_dir: Output directory
        assembly_name: Name of the assembly
        round_num: Racon round number (1-4)
        threads: Number of CPU threads
        state_dir: Directory for state management
        force: Force re-run even if already completed
    
    Returns:
        Path to polished assembly or None if failed
    """
    step = f"racon_{round_num}"
    
    # Check if this round of Racon has already been completed
    if not force and state_dir and is_step_completed(state_dir, assembly_name, step):
        step_data = get_completed_step_data(state_dir, assembly_name, step)
        if step_data and "path" in step_data:
            polished_path = step_data["path"]
            if os.path.exists(polished_path) and os.path.getsize(polished_path) > 0:
                logger.info(f"✅ Resuming: Racon round {round_num} already completed for {assembly_name}")
                return polished_path
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Output file paths
    paf_file = os.path.join(output_dir, f"{assembly_name}.racon{round_num}.paf")
    polished_out = os.path.join(output_dir, f"{assembly_name}.racon{round_num}.fasta")
    
    # Step 1: Map reads to assembly using minimap2
    if os.path.exists(paf_file) and os.path.getsize(paf_file) > 0 and not force:
        logger.info(f"✅ Resuming: Using existing minimap2 alignment for Racon round {round_num}")
    else:
        # Determine minimap2 preset based on read type
        assembly_info = get_assembly_info(state_dir, assembly_name) if state_dir else None
        tech = assembly_info.get("tech") if assembly_info else None
        
        # Set minimap2 preset based on technology type
        preset = "map-pb"  # Default for PacBio
        if tech and tech.startswith("nano"):
            preset = "map-ont"  # For Nanopore
        
        cmd = f"minimap2 -x {preset} -t {threads} {assembly_fasta} {reads_file} > {paf_file}"
        ok, _ = run_cmd(cmd, desc=f"Minimap2 alignment for Racon round {round_num}")
        if not ok or not os.path.exists(paf_file) or os.path.getsize(paf_file) == 0:
            logger.error(f"Minimap2 alignment failed for Racon round {round_num}")
            if state_dir:
                update_state(state_dir, assembly_name, step, "failed", error="minimap2_failed")
            return None
    
    # Step 2: Run Racon for consensus correction
    if os.path.exists(polished_out) and os.path.getsize(polished_out) > 0 and not force:
        logger.info(f"✅ Resuming: Found existing Racon round {round_num} output")
    else:
        cmd = f"racon -t {threads} {reads_file} {paf_file} {assembly_fasta} > {polished_out}"
        ok, _ = run_cmd(cmd, desc=f"Racon polishing round {round_num}")
        if not ok or not os.path.exists(polished_out) or os.path.getsize(polished_out) == 0:
            logger.error(f"Racon polishing round {round_num} failed")
            if state_dir:
                update_state(state_dir, assembly_name, step, "failed", error="racon_failed")
            return None
    
    # Update state
    if state_dir:
        update_state(
            state_dir, assembly_name, step, "completed", 
            path=polished_out, 
            hash=get_file_hash(polished_out),
            timestamp=datetime.now().isoformat()
        )
    
    return polished_out

def run_racon_rounds(assembly_fasta, reads_file, output_dir, assembly_name, rounds=4, threads=16, state_dir=None, force=False, run_stats=True, busco_lineage=None):
    """
    Run multiple rounds of Racon polishing.
    
    Args:
        assembly_fasta: Input assembly FASTA
        reads_file: Long reads file (FASTQ)
        output_dir: Output directory
        assembly_name: Name of the assembly
        rounds: Number of Racon rounds (default: 4)
        threads: Number of CPU threads
        state_dir: Directory for state management
        force: Force re-run even if already completed
        run_stats: Run sequence statistics after polishing
        busco_lineage: BUSCO lineage to use for evaluation (optional)
    
    Returns:
        Path to final polished assembly or None if failed
    """
    # Create polishing directory
    polish_dir = os.path.join(output_dir, f"{assembly_name}_racon")
    os.makedirs(polish_dir, exist_ok=True)
    
    # Store BUSCO results if specified
    if busco_lineage and run_stats:
        busco_dir = os.path.join(output_dir, "busco")
        os.makedirs(busco_dir, exist_ok=True)
    
    # Start with original assembly
    current_assembly = assembly_fasta
    
    # Update original assembly info in state
    if state_dir:
        update_assembly_info(state_dir, assembly_name, original_file=assembly_fasta)
        
        # Get technology type from assembly name or file path
        tech = None
        if "nano" in assembly_name.lower() or "ont" in assembly_name.lower():
            tech = "nano-raw"
        elif "pacbio" in assembly_name.lower() or "pb" in assembly_name.lower():
            if "hifi" in assembly_name.lower():
                tech = "pacbio-hifi"
            else:
                tech = "pacbio-raw"
        
        if tech:
            update_assembly_info(state_dir, assembly_name, tech=tech)
    
    # Run basic stats on original assembly
    if run_stats:
        if state_dir and not is_step_completed(state_dir, assembly_name, "original_stats"):
            stats = run_seqkit_stats(assembly_fasta)
            if stats:
                update_state(state_dir, assembly_name, "original_stats", "completed", result=stats)
                logger.info(f"Original assembly stats for {assembly_name}: N50={stats.get('N50', 'N/A')}, GC={stats.get('gc_content', 'N/A')}")
        
        # Run BUSCO on original assembly if lineage is provided
        if busco_lineage and state_dir and not is_step_completed(state_dir, assembly_name, "original_busco"):
            original_busco_dir = os.path.join(busco_dir, "original")
            run_busco(assembly_fasta, original_busco_dir, busco_lineage, threads, assembly_name, "original", state_dir)
    
    # Run Racon rounds
    for round_num in range(1, rounds + 1):
        logger.info(f"\n=== Racon Round {round_num}/{rounds} for {assembly_name} ===")
        
        # Run this round of Racon
        polished = run_racon_polishing(
            assembly_fasta=current_assembly,
            reads_file=reads_file,
            output_dir=polish_dir,
            assembly_name=assembly_name,
            round_num=round_num,
            threads=threads,
            state_dir=state_dir,
            force=force
        )
        
        if not polished:
            logger.error(f"Racon round {round_num} failed for {assembly_name}")
            return None
        
        # Update current assembly for next round
        current_assembly = polished
        
        # Run stats on final Racon round
        if round_num == rounds and run_stats:
            stats_key = f"racon_{round_num}_stats"
            if state_dir and not is_step_completed(state_dir, assembly_name, stats_key):
                stats = run_seqkit_stats(polished)
                if stats:
                    update_state(state_dir, assembly_name, stats_key, "completed", result=stats)
                    logger.info(f"Racon round {round_num} stats for {assembly_name}: N50={stats.get('N50', 'N/A')}, GC={stats.get('gc_content', 'N/A')}")
            
            # Run BUSCO on final Racon assembly if lineage is provided
            if busco_lineage and state_dir:
                busco_key = f"racon_{round_num}_busco"
                if not is_step_completed(state_dir, assembly_name, busco_key):
                    racon_busco_dir = os.path.join(busco_dir, f"racon_{round_num}")
                    run_busco(polished, racon_busco_dir, busco_lineage, threads, assembly_name, f"racon_{round_num}", state_dir)
    
    return current_assembly

# --------------------------------------------------------------------
# 5. Medaka Polishing (for ONT data)
# --------------------------------------------------------------------
def run_medaka_polishing(assembly_fasta, reads_file, output_dir, assembly_name, model=None, threads=16, state_dir=None, force=False):
    """
    Run Medaka polishing on Nanopore assemblies.
    
    Args:
        assembly_fasta: Input assembly FASTA (typically after Racon)
        reads_file: ONT reads file (FASTQ)
        output_dir: Output directory
        assembly_name: Name of the assembly
        model: Medaka model to use (default: auto-detect)
        threads: Number of CPU threads
        state_dir: Directory for state management
        force: Force re-run even if already completed
    
    Returns:
        Path to Medaka-polished assembly or None if failed
    """
    step = "medaka"
    
    # Check if Medaka has already been completed
    if not force and state_dir and is_step_completed(state_dir, assembly_name, step):
        step_data = get_completed_step_data(state_dir, assembly_name, step)
        if step_data and "path" in step_data:
            polished_path = step_data["path"]
            if os.path.exists(polished_path) and os.path.getsize(polished_path) > 0:
                logger.info(f"✅ Resuming: Medaka polishing already completed for {assembly_name}")
                return polished_path
    
    # Create Medaka output directory
    medaka_dir = os.path.join(output_dir, f"{assembly_name}_medaka")
    os.makedirs(medaka_dir, exist_ok=True)
    
    # Determine the appropriate Medaka model if not specified
    if not model:
        # Default models for different data types
        # Use r941_min_hac_g507 for MinION R9.4.1 high-accuracy basecalling
        model = "r941_min_hac_g507"
        logger.info(f"Auto-selected Medaka model: {model}")
    
    # Output paths
    polished_out = os.path.join(medaka_dir, "consensus.fasta")
    final_out = os.path.join(output_dir, f"{assembly_name}.medaka.fasta")
    
    # Check if Medaka output already exists
    if os.path.exists(polished_out) and os.path.getsize(polished_out) > 0 and not force:
        logger.info(f"✅ Resuming: Found existing Medaka output")
    else:
        # Run Medaka
        cmd = f"medaka_consensus -i {reads_file} -d {assembly_fasta} -o {medaka_dir} -t {threads} -m {model}"
        ok, _ = run_cmd(cmd, desc=f"Medaka polishing")
        if not ok or not os.path.exists(polished_out) or os.path.getsize(polished_out) == 0:
            logger.error(f"Medaka polishing failed for {assembly_name}")
            if state_dir:
                update_state(state_dir, assembly_name, step, "failed", error="medaka_failed")
            return None
    
    # Copy to final output location
    shutil.copy2(polished_out, final_out)
    
    # Update state
    if state_dir:
        update_state(
            state_dir, assembly_name, step, "completed", 
            path=final_out, 
            hash=get_file_hash(final_out),
            timestamp=datetime.now().isoformat()
        )
    
    return final_out

# --------------------------------------------------------------------
# 6. NextPolish (with Illumina short reads)
# --------------------------------------------------------------------
def run_nextpolish(assembly_fasta, r1_file, r2_file, output_dir, assembly_name, threads=16, state_dir=None, force=False):
    """
    Run NextPolish to polish assembly with Illumina short reads.
    
    Args:
        assembly_fasta: Input assembly FASTA (typically after Racon/Medaka)
        r1_file: Illumina R1 reads (FASTQ)
        r2_file: Illumina R2 reads (FASTQ)
        output_dir: Output directory
        assembly_name: Name of the assembly
        threads: Number of CPU threads
        state_dir: Directory for state management
        force: Force re-run even if already completed
    
    Returns:
        Path to NextPolish-polished assembly or None if failed
    """
    step = "nextpolish"
    
    # Check if NextPolish has already been completed
    if not force and state_dir and is_step_completed(state_dir, assembly_name, step):
        step_data = get_completed_step_data(state_dir, assembly_name, step)
        if step_data and "path" in step_data:
            polished_path = step_data["path"]
            if os.path.exists(polished_path) and os.path.getsize(polished_path) > 0:
                logger.info(f"✅ Resuming: NextPolish already completed for {assembly_name}")
                return polished_path
    
    # Create NextPolish output directory
    nextpolish_dir = os.path.join(output_dir, f"{assembly_name}_nextpolish")
    os.makedirs(nextpolish_dir, exist_ok=True)
    
    # Prepare NextPolish input files
    input_genome = os.path.join(nextpolish_dir, "input.fasta")
    shutil.copy2(assembly_fasta, input_genome)
    
    # Create a sgs.fofn file with the Illumina read paths
    sgs_fofn = os.path.join(nextpolish_dir, "sgs.fofn")
    with open(sgs_fofn, 'w') as f:
        f.write(f"{r1_file}\n{r2_file}\n")
    
    # Create NextPolish configuration file
    config_file = os.path.join(nextpolish_dir, "run.cfg")
    with open(config_file, 'w') as f:
        f.write(f"""
[General]
job_type = local
job_prefix = nextpolish
task = best
rewrite = yes
rerun = 3
parallel_jobs = {threads}
multithread_jobs = {threads}
genome = {input_genome}
genome_size = auto
workdir = {nextpolish_dir}
polish_options = -p {threads}

[sgs_option]
sgs_fofn = {sgs_fofn}
sgs_options = -max_depth 100 -bwa
        """)
    
    # Output paths
    polished_out = os.path.join(nextpolish_dir, "genome.nextpolish.fasta")
    final_out = os.path.join(output_dir, f"{assembly_name}.nextpolish.fasta")
    
    # Check if NextPolish output already exists
    if os.path.exists(polished_out) and os.path.getsize(polished_out) > 0 and not force:
        logger.info(f"✅ Resuming: Found existing NextPolish output")
    else:
        # Run NextPolish
        cmd = f"cd {nextpolish_dir} && nextPolish {config_file}"
        ok, _ = run_cmd(cmd, desc=f"NextPolish with Illumina reads", cwd=nextpolish_dir)
        if not ok or not os.path.exists(polished_out) or os.path.getsize(polished_out) == 0:
            logger.error(f"NextPolish failed for {assembly_name}")
            # Look for alternative output location
            alt_out = os.path.join(nextpolish_dir, "genome.nextpolish")
            if os.path.exists(alt_out) and os.path.getsize(alt_out) > 0:
                polished_out = alt_out
            else:
                if state_dir:
                    update_state(state_dir, assembly_name, step, "failed", error="nextpolish_failed")
                return None
    
    # Copy to final output location
    shutil.copy2(polished_out, final_out)
    
    # Update state
    if state_dir:
        update_state(
            state_dir, assembly_name, step, "completed", 
            path=final_out, 
            hash=get_file_hash(final_out),
            timestamp=datetime.now().isoformat()
        )
    
    return final_out

# --------------------------------------------------------------------
# 7. Complete Polishing Pipeline
# --------------------------------------------------------------------
def run_full_polishing_pipeline(assembly_fasta, long_reads, r1_file, r2_file, output_dir, 
                               assembly_name, tech_type=None, racon_rounds=4, threads=16, 
                               medaka_model=None, busco_lineage=None, force=False):
    """
    Run the complete polishing pipeline: Racon -> Medaka (ONT only) -> NextPolish.
    
    Args:
        assembly_fasta: Input assembly FASTA
        long_reads: Long reads file (FASTQ) for Racon
        r1_file: Illumina R1 reads for NextPolish
        r2_file: Illumina R2 reads for NextPolish
        output_dir: Output directory
        assembly_name: Name of the assembly
        tech_type: Technology type (pacbio-hifi, pacbio-raw, nano-raw, etc.)
        racon_rounds: Number of Racon rounds
        threads: Number of CPU threads
        medaka_model: Medaka model for ONT data
        busco_lineage: BUSCO lineage for evaluation
        force: Force re-run all steps
    
    Returns:
        Dictionary with paths to polished assemblies at each stage
    """
    logger.info(f"\n=== Starting Polishing Pipeline for {assembly_name} ===")
    logger.info(f"Input assembly: {assembly_fasta}")
    logger.info(f"Long reads: {long_reads}")
    logger.info(f"Illumina R1: {r1_file}")
    logger.info(f"Illumina R2: {r2_file}")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Update assembly info in state
    if tech_type:
        update_assembly_info(output_dir, assembly_name, tech=tech_type, original_file=assembly_fasta)
    
    # Results dictionary
    results = {
        "original": assembly_fasta,
        "racon": None,
        "medaka": None,
        "nextpolish": None
    }
    
    # 1. Run Racon rounds
    logger.info(f"\n=== Racon Polishing ({racon_rounds} rounds) ===")
    racon_out = run_racon_rounds(
        assembly_fasta=assembly_fasta,
        reads_file=long_reads,
        output_dir=output_dir,
        assembly_name=assembly_name,
        rounds=racon_rounds,
        threads=threads,
        state_dir=output_dir,
        force=force,
        run_stats=True,
        busco_lineage=busco_lineage
    )
    
    if not racon_out:
        logger.error("Racon polishing failed. Stopping pipeline.")
        return results
    
    results["racon"] = racon_out
    current_assembly = racon_out
    
    # 2. Run Medaka (for ONT data only)
    is_ont = tech_type and ("nano" in tech_type or "ont" in tech_type.lower())
    is_ont = is_ont or ("nano" in assembly_name.lower() or "ont" in assembly_name.lower())
    
    if is_ont:
        logger.info(f"\n=== Medaka Polishing (ONT data) ===")
        medaka_out = run_medaka_polishing(
            assembly_fasta=current_assembly,
            reads_file=long_reads,
            output_dir=output_dir,
            assembly_name=assembly_name,
            model=medaka_model,
            threads=threads,
            state_dir=output_dir,
            force=force
        )
        
        if medaka_out:
            results["medaka"] = medaka_out
            current_assembly = medaka_out
            
            # Run stats on Medaka output
            if not is_step_completed(output_dir, assembly_name, "medaka_stats"):
                stats = run_seqkit_stats(medaka_out)
                if stats:
                    update_state(output_dir, assembly_name, "medaka_stats", "completed", result=stats)
                    logger.info(f"Medaka stats for {assembly_name}: N50={stats.get('N50', 'N/A')}, GC={stats.get('gc_content', 'N/A')}")
            
            # Run BUSCO on Medaka output if lineage is provided
            if busco_lineage and not is_step_completed(output_dir, assembly_name, "medaka_busco"):
                busco_dir = os.path.join(output_dir, "busco", "medaka")
                run_busco(medaka_out, busco_dir, busco_lineage, threads, assembly_name, "medaka", output_dir)
        else:
            logger.warning("Medaka polishing failed. Continuing with Racon output.")
    else:
        logger.info("Skipping Medaka polishing (not ONT data)")
    
    # 3. Run NextPolish with Illumina reads
    if r1_file and r2_file and os.path.exists(r1_file) and os.path.exists(r2_file):
        logger.info(f"\n=== NextPolish with Illumina Reads ===")
        nextpolish_out = run_nextpolish(
            assembly_fasta=current_assembly,
            r1_file=r1_file,
            r2_file=r2_file,
            output_dir=output_dir,
            assembly_name=assembly_name,
            threads=threads,
            state_dir=output_dir,
            force=force
        )
        
        if nextpolish_out:
            results["nextpolish"] = nextpolish_out
            
            # Run stats on NextPolish output
            if not is_step_completed(output_dir, assembly_name, "nextpolish_stats"):
                stats = run_seqkit_stats(nextpolish_out)
                if stats:
                    update_state(output_dir, assembly_name, "nextpolish_stats", "completed", result=stats)
                    logger.info(f"NextPolish stats for {assembly_name}: N50={stats.get('N50', 'N/A')}, GC={stats.get('gc_content', 'N/A')}")
            
            # Run BUSCO on NextPolish output if lineage is provided
            if busco_lineage and not is_step_completed(output_dir, assembly_name, "nextpolish_busco"):
                busco_dir = os.path.join(output_dir, "busco", "nextpolish")
                run_busco(nextpolish_out, busco_dir, busco_lineage, threads, assembly_name, "nextpolish", output_dir)
        else:
            logger.warning("NextPolish failed. Final output will be from previous step.")
    else:
        logger.warning("Skipping NextPolish (Illumina reads not provided or not found)")
    
    # Generate a summary of the polishing results
    create_polishing_summary(output_dir, assembly_name)
    
    return results

# --------------------------------------------------------------------
# 8. Summary and Evaluation
# --------------------------------------------------------------------
def create_polishing_summary(output_dir, assembly_name=None):
    """
    Create a summary of polishing results, showing improvements at each stage.
    """
    state = read_state(output_dir)
    results = []
    
    # Filter for specific assembly if provided
    if assembly_name:
        if assembly_name in state:
            assemblies = {assembly_name: state[assembly_name]}
        else:
            logger.warning(f"No data found for assembly {assembly_name}")
            return
    else:
        assemblies = state
    
    # Collect results for each assembly
    for name, data in assemblies.items():
        # Get stats for each stage
        stages = ["original", "racon_4", "medaka", "nextpolish"]
        row = {"assembly": name}
        
        # Get path to each polished assembly
        polishing_steps = data.get("polishing_steps", {})
        for step_name, step_info in polishing_steps.items():
            if step_info.get("status") == "completed" and "path" in step_info:
                stage = step_name
                row[f"{stage}_path"] = step_info["path"]
        
        # Get stats for each stage
        evaluations = data.get("evaluations", {})
        for stage in stages:
            stats_key = f"{stage}_stats"
            busco_key = f"{stage}_busco"
            
            # Add basic stats if available
            if stats_key in evaluations and "result" in evaluations[stats_key]:
                stats = evaluations[stats_key]["result"]
                row[f"{stage}_total_bp"] = stats.get("sum_len", 0)
                row[f"{stage}_contigs"] = stats.get("num_seqs", 0)
                row[f"{stage}_n50"] = stats.get("N50", 0)
                row[f"{stage}_gc"] = stats.get("gc_content", 0.0)
            
            # Add BUSCO results if available
            if busco_key in evaluations and "result" in evaluations[busco_key]:
                busco = evaluations[busco_key]["result"]
                row[f"{stage}_busco"] = busco.get("complete_percent", 0.0)
                row[f"{stage}_missing"] = busco.get("missing_percent", 0.0)
        
        results.append(row)
    
    if not results:
        logger.warning("No polishing results found")
        return
    
    # Create a DataFrame
    df = pd.DataFrame(results)
    
    # Calculate improvements
    if len(results) > 0 and "original_busco" in df.columns and "nextpolish_busco" in df.columns:
        df["busco_improvement"] = df["nextpolish_busco"] - df["original_busco"]
    if len(results) > 0 and "original_n50" in df.columns and "nextpolish_n50" in df.columns:
        df["n50_improvement"] = df["nextpolish_n50"] - df["original_n50"]
    
    # Save to CSV
    csv_path = os.path.join(output_dir, "polishing_summary.csv")
    df.to_csv(csv_path, index=False)
    
    # Save markdown summary
    md_path = os.path.join(output_dir, "polishing_summary.md")
    with open(md_path, 'w') as f:
        f.write("# Assembly Polishing Summary\n\n")
        f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        # For each assembly, create a section with detailed results
        for name, data in assemblies.items():
            f.write(f"## {name}\n\n")
            
            # Create a table with key metrics
            f.write("| Stage | Contigs | Total bp | N50 | GC% | BUSCO Complete% | BUSCO Missing% |\n")
            f.write("|-------|---------|----------|-----|-----|-----------------|----------------|\n")
            
            # Get stats for each stage
            stages = ["original", "racon_4", "medaka", "nextpolish"]
            evaluations = data.get("evaluations", {})
            
            for stage in stages:
                stats_key = f"{stage}_stats"
                busco_key = f"{stage}_busco"
                
                # Get stats and BUSCO results
                stats = evaluations.get(stats_key, {}).get("result", {})
                busco = evaluations.get(busco_key, {}).get("result", {})
                
                # Format row
                contigs = stats.get("num_seqs", "N/A")
                total_bp = f"{stats.get('sum_len', 0):,}"
                n50 = f"{stats.get('N50', 0):,}"
                gc = f"{stats.get('gc_content', 0):.2f}%"
                busco_complete = f"{busco.get('complete_percent', 0):.2f}%"
                busco_missing = f"{busco.get('missing_percent', 0):.2f}%"
                
                f.write(f"| {stage.capitalize()} | {contigs} | {total_bp} | {n50} | {gc} | {busco_complete} | {busco_missing} |\n")
            
            f.write("\n")
            
            # Create a summary of improvements
            if "nextpolish_stats" in evaluations and "original_stats" in evaluations:
                original_stats = evaluations["original_stats"].get("result", {})
                final_stats = evaluations["nextpolish_stats"].get("result", {})
                
                if original_stats and final_stats:
                    n50_change = final_stats.get("N50", 0) - original_stats.get("N50", 0)
                    n50_pct = (n50_change / original_stats.get("N50", 1)) * 100 if original_stats.get("N50", 0) > 0 else 0
                    
                    f.write("### N50 Improvement\n")
                    f.write(f"- Original N50: {original_stats.get('N50', 0):,} bp\n")
                    f.write(f"- Final N50: {final_stats.get('N50', 0):,} bp\n")
                    f.write(f"- Change: {n50_change:+,} bp ({n50_pct:+.2f}%)\n\n")
            
            if "nextpolish_busco" in evaluations and "original_busco" in evaluations:
                original_busco = evaluations["original_busco"].get("result", {})
                final_busco = evaluations["nextpolish_busco"].get("result", {})
                
                if original_busco and final_busco:
                    busco_change = final_busco.get("complete_percent", 0) - original_busco.get("complete_percent", 0)
                    
                    f.write("### BUSCO Improvement\n")
                    f.write(f"- Original BUSCO complete: {original_busco.get('complete_percent', 0):.2f}%\n")
                    f.write(f"- Final BUSCO complete: {final_busco.get('complete_percent', 0):.2f}%\n")
                    f.write(f"- Change: {busco_change:+.2f}%\n\n")
            
            # Add paths to final assemblies
            f.write("### Assembly Paths\n")
            polishing_steps = data.get("polishing_steps", {})
            for step_name, step_info in polishing_steps.items():
                if step_info.get("status") == "completed" and "path" in step_info:
                    f.write(f"- {step_name.capitalize()}: `{step_info['path']}`\n")
            
            f.write("\n---\n\n")
    
    logger.info(f"Polishing summary saved to:\n  {csv_path}\n  {md_path}")

def check_dependencies():
    """
    Check for required command-line tools.
    Exits if any are missing.
    """
    import shutil
    
    required_tools = [
        "minimap2",
        "racon",
        "medaka_consensus",
        "nextPolish",
        "seqkit",
        "busco"
    ]
    
    missing_tools = [tool for tool in required_tools if shutil.which(tool) is None]
    
    if missing_tools:
        print("\n❌ Missing dependencies:")
        print("🔧 Missing tools:")
        for tool in missing_tools:
            print(f"  - {tool}")
        
        print("\n💡 Please install the missing components and try again.")
        print("💡 If you're using Conda/Mamba, you can install most tools with:")
        print("    conda install -c bioconda minimap2 racon medaka nextpolish busco seqkit")
        print()
        sys.exit(1)
    
    print("✅ All dependencies satisfied.\n")

# --------------------------------------------------------------------
# 9. Main
# --------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Assembly Polishing Pipeline (Racon + Medaka + NextPolish)"
    )
    
    # Required inputs
    parser.add_argument("-a", "--assembly", required=True, help="Input assembly FASTA file")
    parser.add_argument("-o", "--output_dir", default="Polished_Assemblies", help="Output directory")
    parser.add_argument("-l", "--long_reads", required=True, help="Long reads FASTQ/FASTA file for Racon/Medaka")
    parser.add_argument("-1", "--r1", help="Illumina R1 reads for NextPolish")
    parser.add_argument("-2", "--r2", help="Illumina R2 reads for NextPolish")
    
    # Optional parameters
    parser.add_argument("-t", "--tech", default=None, choices=["pacbio-hifi", "pacbio-raw", "nano-raw", "nano-hq", "nano-corr"], 
                        help="Long read technology type")
    parser.add_argument("-n", "--name", help="Assembly name (defaults to filename)")
    parser.add_argument("--threads", type=int, default=16, help="Number of CPU threads")
    parser.add_argument("--racon_rounds", type=int, default=4, help="Number of Racon rounds (default: 4)")
    parser.add_argument("--medaka_model", default=None, help="Medaka model (default: auto-detect)")
    parser.add_argument("--busco", default=None, help="BUSCO lineage for evaluation (e.g., mollusca_odb10)")
    
    # Control flags
    parser.add_argument("--force", action="store_true", help="Force re-run all steps")
    parser.add_argument("--force-racon", action="store_true", help="Force re-run Racon steps")
    parser.add_argument("--force-medaka", action="store_true", help="Force re-run Medaka step")
    parser.add_argument("--force-nextpolish", action="store_true", help="Force re-run NextPolish step")
    parser.add_argument("--skip-medaka", action="store_true", help="Skip Medaka polishing even for ONT data")
    parser.add_argument("--skip-nextpolish", action="store_true", help="Skip NextPolish step")
    parser.add_argument("--show-state", action="store_true", help="Show current state information and exit")
    parser.add_argument("--reset-state", action="store_true", help="Reset pipeline state before running")
    
    args = parser.parse_args()
    
    # Setup logging
    log_fp = setup_logging(args.output_dir)
    logger.info(f"=== Assembly Polishing Pipeline: Racon + Medaka + NextPolish ===")
    logger.info(f"Input assembly: {args.assembly}")
    logger.info(f"Long reads: {args.long_reads}")
    logger.info(f"Illumina R1: {args.r1 if args.r1 else 'Not provided'}")
    logger.info(f"Illumina R2: {args.r2 if args.r2 else 'Not provided'}")
    logger.info(f"Output dir: {args.output_dir}")
    logger.info(f"Threads: {args.threads}")
    logger.info(f"Logs -> {log_fp}")
    
    # Check dependencies
    check_dependencies()
    
    # Determine assembly name if not provided
    if not args.name:
        args.name = os.path.basename(args.assembly).split('.')[0]
    
    # Reset state if requested
    if args.reset_state:
        state_path = state_file_path(args.output_dir)
        if os.path.exists(state_path):
            os.unlink(state_path)
            logger.info(f"🧹 Reset pipeline state file: {state_path}")
    
    # Show state if requested
    if args.show_state:
        state = read_state(args.output_dir)
        if not state:
            logger.info("No pipeline state found or state is empty.")
        else:
            logger.info("Current pipeline state:")
            for assembly, data in state.items():
                logger.info(f"Assembly: {assembly}")
                
                # Show polishing steps
                if "polishing_steps" in data:
                    logger.info("  Polishing steps:")
                    for step, step_data in data["polishing_steps"].items():
                        status = step_data.get("status", "unknown")
                        status_symbol = "✅" if status == "completed" else "❌"
                        logger.info(f"    {step}: {status_symbol} {status}")
                
                # Show evaluations
                if "evaluations" in data:
                    logger.info("  Evaluations:")
                    for eval_step, eval_data in data["evaluations"].items():
                        status = eval_data.get("status", "unknown")
                        status_symbol = "✅" if status == "completed" else "❌"
                        logger.info(f"    {eval_step}: {status_symbol} {status}")
            
            # Show summary
            create_polishing_summary(args.output_dir)
        
        sys.exit(0)
    
    # Create necessary directories
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Update assembly information in state
    update_assembly_info(args.output_dir, args.name, tech=args.tech, original_file=args.assembly)
    
    # Determine force flags for each step
    force_racon = args.force or args.force_racon
    force_medaka = args.force or args.force_medaka
    force_nextpolish = args.force or args.force_nextpolish
    
    # Run Racon rounds
    current_assembly = args.assembly
    if not args.show_state:
        logger.info(f"\n=== Racon Polishing ({args.racon_rounds} rounds) ===")
        racon_out = run_racon_rounds(
            assembly_fasta=args.assembly,
            reads_file=args.long_reads,
            output_dir=args.output_dir,
            assembly_name=args.name,
            rounds=args.racon_rounds,
            threads=args.threads,
            state_dir=args.output_dir,
            force=force_racon,
            run_stats=True,
            busco_lineage=args.busco
        )
        
        if racon_out:
            current_assembly = racon_out
            logger.info(f"Racon polishing completed: {racon_out}")
        else:
            logger.error("Racon polishing failed")
            sys.exit(1)
    
    # Run Medaka (for ONT data only)
    if not args.skip_medaka:
        # Check if this is ONT data
        is_ont = args.tech and ("nano" in args.tech or "ont" in args.tech.lower())
        is_ont = is_ont or ("nano" in args.name.lower() or "ont" in args.name.lower())
        
        if is_ont:
            logger.info(f"\n=== Medaka Polishing (ONT data) ===")
            medaka_out = run_medaka_polishing(
                assembly_fasta=current_assembly,
                reads_file=args.long_reads,
                output_dir=args.output_dir,
                assembly_name=args.name,
                model=args.medaka_model,
                threads=args.threads,
                state_dir=args.output_dir,
                force=force_medaka
            )
            
            if medaka_out:
                current_assembly = medaka_out
                logger.info(f"Medaka polishing completed: {medaka_out}")
                
                # Run stats and BUSCO on Medaka output
                if not is_step_completed(args.output_dir, args.name, "medaka_stats"):
                    stats = run_seqkit_stats(medaka_out)
                    if stats:
                        update_state(args.output_dir, args.name, "medaka_stats", "completed", result=stats)
                
                if args.busco and not is_step_completed(args.output_dir, args.name, "medaka_busco"):
                    busco_dir = os.path.join(args.output_dir, "busco", "medaka")
                    run_busco(medaka_out, busco_dir, args.busco, args.threads, args.name, "medaka", args.output_dir)
            else:
                logger.warning("Medaka polishing failed. Continuing with Racon output.")
        else:
            logger.info("Skipping Medaka polishing (not ONT data)")
    
    # Run NextPolish with Illumina reads
    if not args.skip_nextpolish and args.r1 and args.r2 and os.path.exists(args.r1) and os.path.exists(args.r2):
        logger.info(f"\n=== NextPolish with Illumina Reads ===")
        nextpolish_out = run_nextpolish(
            assembly_fasta=current_assembly,
            r1_file=args.r1,
            r2_file=args.r2,
            output_dir=args.output_dir,
            assembly_name=args.name,
            threads=args.threads,
            state_dir=args.output_dir,
            force=force_nextpolish
        )
        
        if nextpolish_out:
            current_assembly = nextpolish_out
            logger.info(f"NextPolish completed: {nextpolish_out}")
            
            # Run stats and BUSCO on NextPolish output
            if not is_step_completed(args.output_dir, args.name, "nextpolish_stats"):
                stats = run_seqkit_stats(nextpolish_out)
                if stats:
                    update_state(args.output_dir, args.name, "nextpolish_stats", "completed", result=stats)
            
            if args.busco and not is_step_completed(args.output_dir, args.name, "nextpolish_busco"):
                busco_dir = os.path.join(args.output_dir, "busco", "nextpolish")
                run_busco(nextpolish_out, busco_dir, args.busco, args.threads, args.name, "nextpolish", args.output_dir)
        else:
            logger.warning("NextPolish failed. Final output will be from previous step.")
    else:
        logger.info("Skipping NextPolish (Illumina reads not provided or not found)")
    
    # Generate a summary of the polishing results
    create_polishing_summary(args.output_dir, args.name)
    
    logger.info(f"\n=== Polishing Pipeline Completed ===")
    logger.info(f"Final polished assembly: {current_assembly}")
    logger.info(f"Summary reports created in {args.output_dir}")

if __name__ == "__main__":
    main()