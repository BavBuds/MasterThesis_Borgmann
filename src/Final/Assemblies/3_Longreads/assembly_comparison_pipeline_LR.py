#!/usr/bin/env python3
"""
Long-read Assembly Pipeline (FLYE vs CANU)
with BUSCO (mollusca_odb10 & eukaryota_odb10) and GC% via seqkit.
Includes robust resume functionality to continue from interrupted runs.

Author: Based on Max Borgmann's pipeline
Date: 03.04.2025

Steps:
 1. For each sample:
    A) Run FLYE for PacBio HiFi
    B) Run FLYE for ONT
    C) Run CANU for PacBio HiFi
    D) Run CANU for ONT
 2. After each assembly, run BUSCO *twice*:
    - with mollusca_odb10
    - with eukaryota_odb10
    Parse the JSON to get:
       - complete_percent
       - Contigs N50
       - # contigs, total length, missing%
 3. Also run seqkit stats to get GC%.
 4. Summarize everything in a single table with columns:
    - busco_mollusca
    - busco_eukaryota
    - n50
    - missing_mollusca
    - missing_eukaryota
    - num_contigs
    - total_length
    - gc
    Then sort by [busco_mollusca DESC, busco_eukaryota DESC, n50 DESC].
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
logger = logging.getLogger("longread_assembly_comparison")
logger.setLevel(logging.DEBUG)

def setup_logging(log_dir, log_name="longread_assembly_comparison.log"):
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
    return os.path.join(output_dir, "pipeline_state.json")

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
        "sample_name": {
            "last_update": "timestamp",
            "assemblers": {
                "FLYE_pacbio-hifi": {
                    "assembly": {"status": "completed", "path": "...", "hash": "..."},
                    "busco_mollusca": {"status": "completed", "path": "...", "result": {...}},
                    "busco_eukaryota": {"status": "completed", "path": "...", "result": {...}},
                    "gc": {"status": "completed", "value": 0.42}
                },
                # Other assemblers...
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

def update_state(output_dir, sample, assembler, config, step, status, **kwargs):
    """
    Update the pipeline state for a specific step.
    
    Args:
        output_dir: Output directory containing the state file
        sample: Sample name
        assembler: Assembler name (FLYE, CANU)
        config: Configuration (pacbio-hifi, nano-raw, etc.)
        step: Step name (assembly, busco_mollusca, busco_eukaryota, gc)
        status: Status of the step (completed, failed)
        **kwargs: Additional data to store (path, hash, result, value, etc.)
    """
    state = read_state(output_dir)
    
    # Ensure nested structure exists
    if sample not in state:
        state[sample] = {"last_update": datetime.now().isoformat(), "assemblers": {}}
    
    assembler_key = f"{assembler}_{config}"
    if assembler_key not in state[sample]["assemblers"]:
        state[sample]["assemblers"][assembler_key] = {}
    
    # Update the step
    state[sample]["assemblers"][assembler_key][step] = {"status": status, **kwargs}
    state[sample]["last_update"] = datetime.now().isoformat()
    
    # Write the updated state
    write_state(state, output_dir)

def is_step_completed(output_dir, sample, assembler, config, step, verify_file=None):
    """
    Check if a step has been completed successfully.
    If verify_file is provided, also check if the file exists and hash matches.
    """
    state = read_state(output_dir)
    
    # Check if step exists in state and was completed
    try:
        step_state = state[sample]["assemblers"][f"{assembler}_{config}"][step]
        if step_state["status"] != "completed":
            return False
        
        # If a file path is given for verification
        if verify_file and "path" in step_state:
            stored_path = step_state["path"]
            
            # Check if the file exists
            if not os.path.exists(stored_path):
                logger.warning(f"File marked as completed but doesn't exist: {stored_path}")
                return False
            
            # If hash is stored, verify the hash
            if "hash" in step_state:
                current_hash = get_file_hash(stored_path)
                if current_hash != step_state["hash"]:
                    logger.warning(f"File hash mismatch for {stored_path}")
                    return False
        
        return True
    except (KeyError, TypeError):
        return False

def get_completed_step_data(output_dir, sample, assembler, config, step):
    """Retrieve data for a completed step from the state file."""
    state = read_state(output_dir)
    try:
        return state[sample]["assemblers"][f"{assembler}_{config}"][step]
    except (KeyError, TypeError):
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

def get_gc_content(fasta_file, out_file=None, output_dir=None, sample=None, assembler=None, config=None, force=False):
    """
    Use seqkit to get GC% of the entire assembly.
    Implements resume functionality.
    """
    if not fasta_file or not os.path.exists(fasta_file):
        return 0.0
    
    # Check if we've already calculated GC content
    if output_dir and sample and assembler and config and not force:
        step_data = get_completed_step_data(output_dir, sample, assembler, config, "gc")
        if step_data and "value" in step_data:
            logger.info(f"✅ Resuming: GC content already calculated for {sample}/{assembler}/{config}")
            return step_data["value"]
    
    if not out_file:
        out_file = f"{fasta_file}.seqkit_stats.txt"
    
    # If stats file already exists and not force, use it
    if os.path.exists(out_file) and not force:
        logger.info(f"✅ Resuming: Using existing seqkit stats file: {out_file}")
    else:
        cmd = f"seqkit stats -a {fasta_file} -o {out_file}"
        ok, _ = run_cmd(cmd, desc="seqkit stats")
        if not ok or not os.path.exists(out_file):
            return 0.0

    df = pd.read_csv(out_file, sep='\t', comment='#')
    if df.empty:
        return 0.0
    
    gc_value = float(df.iloc[0].get("gc_content", 0.0))
    
    # Update state
    if output_dir and sample and assembler and config:
        update_state(output_dir, sample, assembler, config, "gc", "completed", value=gc_value)
    
    return gc_value

# --------------------------------------------------------------------
# 4. FLYE Assembly
# --------------------------------------------------------------------
def run_flye_assembly(reads, output_dir, sample_name, read_type, genome_size, threads=16, 
                      polishing_iterations=1, state_dir=None, force=False):
    """
    Run FLYE assembler with appropriate parameters based on read type.
    Implements resume functionality.
    
    Parameters:
        reads: Path to long read file (FASTQ/FASTA)
        output_dir: Output directory
        sample_name: Sample name prefix
        read_type: One of "pacbio-hifi", "pacbio-raw", "nano-raw", "nano-hq", "nano-corr"
        genome_size: Approximate genome size (e.g., "50m" for 50 Mbp)
        threads: Number of CPU threads
        polishing_iterations: Number of polishing rounds
        state_dir: Directory containing state file for resume functionality
        force: Force re-run even if previously completed
    
    Returns:
        Path to final assembly FASTA or None if failed
    """
    assembler = "FLYE"
    config = read_type
    
    # Check if this assembly was already completed
    if not force and state_dir and is_step_completed(state_dir, sample_name, assembler, config, "assembly"):
        step_data = get_completed_step_data(state_dir, sample_name, assembler, config, "assembly")
        if step_data and "path" in step_data:
            assembly_path = step_data["path"]
            if os.path.exists(assembly_path) and os.path.getsize(assembly_path) > 0:
                logger.info(f"✅ Resuming: FLYE assembly ({read_type}) already completed for {sample_name}")
                return assembly_path
    
    os.makedirs(output_dir, exist_ok=True)
    flye_dir = os.path.join(output_dir, f"{sample_name}_flye_{read_type}")
    os.makedirs(flye_dir, exist_ok=True)
    
    # Map user-friendly types to Flye parameters
    tech_map = {
        "pacbio-hifi": "--pacbio-hifi",
        "pacbio-raw": "--pacbio-raw",
        "nano-raw": "--nano-raw",
        "nano-hq": "--nano-hq",
        "nano-corr": "--nano-corr"
    }
    
    if read_type not in tech_map:
        logger.error(f"Invalid read_type for FLYE: {read_type}")
        logger.error(f"Must be one of: {', '.join(tech_map.keys())}")
        if state_dir:
            update_state(state_dir, sample_name, assembler, config, "assembly", "failed", 
                         error="invalid_read_type")
        return None
    
    tech_param = tech_map[read_type]
    
    # Check if final assembly already exists and we're not forcing a re-run
    final_asm = os.path.join(flye_dir, "assembly.fasta")
    if os.path.exists(final_asm) and os.path.getsize(final_asm) > 0 and not force:
        logger.info(f"✅ Resuming: Found existing FLYE assembly: {final_asm}")
    else:
        # Run FLYE assembly
        cmd = (
            f"flye {tech_param} {reads} --out-dir {flye_dir} --threads {threads} "
            f"--genome-size {genome_size} --iterations {polishing_iterations}"
        )
        ok, _ = run_cmd(cmd, desc=f"FLYE assembly ({read_type})")
        if not ok:
            if state_dir:
                update_state(state_dir, sample_name, assembler, config, "assembly", "failed",
                             error="flye_assembly_failed")
            return None
    
        # Check for final assembly file
        if not os.path.exists(final_asm):
            logger.error(f"FLYE final assembly not found: {final_asm}")
            if state_dir:
                update_state(state_dir, sample_name, assembler, config, "assembly", "failed",
                             error="flye_missing_assembly")
            return None
    
    # Copy the result to the output directory
    out_fa = os.path.join(output_dir, f"{sample_name}_flye_{read_type}.fasta")
    shutil.copy2(final_asm, out_fa)
    
    # Update state
    if state_dir:
        update_state(
            state_dir, sample_name, assembler, config, "assembly", "completed", 
            path=out_fa, 
            hash=get_file_hash(out_fa),
            timestamp=datetime.now().isoformat()
        )
    
    return out_fa

# --------------------------------------------------------------------
# 5. CANU Assembly
# --------------------------------------------------------------------
def run_canu_assembly(reads, output_dir, sample_name, read_type, genome_size, threads=16,
                      state_dir=None, force=False):
    """
    Run CANU assembler with appropriate parameters based on read type.
    Implements resume functionality.
    
    Parameters:
        reads: Path to long read file (FASTQ/FASTA)
        output_dir: Output directory
        sample_name: Sample name prefix
        read_type: One of "pacbio-hifi", "pacbio-raw", "nanopore-raw", "nanopore-corr"
        genome_size: Approximate genome size (e.g., "50m" for 50 Mbp)
        threads: Number of CPU threads
        state_dir: Directory containing state file for resume functionality
        force: Force re-run even if previously completed
    
    Returns:
        Path to final assembly FASTA or None if failed
    """
    assembler = "CANU"
    config = read_type
    
    # Check if this assembly was already completed
    if not force and state_dir and is_step_completed(state_dir, sample_name, assembler, config, "assembly"):
        step_data = get_completed_step_data(state_dir, sample_name, assembler, config, "assembly")
        if step_data and "path" in step_data:
            assembly_path = step_data["path"]
            if os.path.exists(assembly_path) and os.path.getsize(assembly_path) > 0:
                logger.info(f"✅ Resuming: CANU assembly ({read_type}) already completed for {sample_name}")
                return assembly_path
    
    os.makedirs(output_dir, exist_ok=True)
    canu_dir = os.path.join(output_dir, f"{sample_name}_canu_{read_type}")
    os.makedirs(canu_dir, exist_ok=True)
    
    # Map user-friendly types to Canu parameters
    tech_map = {
        "pacbio-hifi": "-pacbio-hifi",
        "pacbio-raw": "-pacbio",
        "nanopore-raw": "-nanopore",
        "nanopore-corr": "-nanopore-corrected"
    }
    
    if read_type not in tech_map:
        logger.error(f"Invalid read_type for CANU: {read_type}")
        logger.error(f"Must be one of: {', '.join(tech_map.keys())}")
        if state_dir:
            update_state(state_dir, sample_name, assembler, config, "assembly", "failed",
                         error="invalid_read_type")
        return None
    
    tech_param = tech_map[read_type]
    
    # Set some reasonable error rates based on read type
    error_rates = {
        "pacbio-hifi": "0.01",
        "pacbio-raw": "0.045",
        "nanopore-raw": "0.15",
        "nanopore-corr": "0.05"
    }
    
    error_rate = error_rates.get(read_type, "0.045")
    
    # Check if final assembly already exists
    final_asm = os.path.join(canu_dir, f"{sample_name}.contigs.fasta")
    if os.path.exists(final_asm) and os.path.getsize(final_asm) > 0 and not force:
        logger.info(f"✅ Resuming: Found existing CANU assembly: {final_asm}")
    else:
        # Run CANU assembly
        prefix = os.path.join(canu_dir, sample_name)
        cmd = (
            f"canu -p {sample_name} -d {canu_dir} {tech_param} {reads} "
            f"genomeSize={genome_size} errorRate={error_rate} "
            f"useGrid=false maxThreads={threads} maxMemory=32g"
        )
        ok, _ = run_cmd(cmd, desc=f"CANU assembly ({read_type})")
        if not ok:
            if state_dir:
                update_state(state_dir, sample_name, assembler, config, "assembly", "failed",
                             error="canu_assembly_failed")
            return None
    
        # Check for final assembly file
        if not os.path.exists(final_asm):
            logger.error(f"CANU final assembly not found: {final_asm}")
            if state_dir:
                update_state(state_dir, sample_name, assembler, config, "assembly", "failed",
                             error="canu_missing_assembly")
            return None
    
    # Copy the result to the output directory
    out_fa = os.path.join(output_dir, f"{sample_name}_canu_{read_type}.fasta")
    shutil.copy2(final_asm, out_fa)
    
    # Update state
    if state_dir:
        update_state(
            state_dir, sample_name, assembler, config, "assembly", "completed", 
            path=out_fa, 
            hash=get_file_hash(out_fa),
            timestamp=datetime.now().isoformat()
        )
    
    return out_fa

# --------------------------------------------------------------------
# 6. BUSCO (run for two lineages)
# --------------------------------------------------------------------
def run_busco(fasta_file, out_dir, lineage, threads, sample, assembler, config, state_dir=None, force=False):
    """
    Runs BUSCO on 'fasta_file' with the given 'lineage' with resume functionality.
    """
    step = f"busco_{lineage.split('_')[0]}"  # e.g., busco_mollusca
    
    # Check if BUSCO has already been run
    if not force and state_dir and is_step_completed(state_dir, sample, assembler, config, step):
        step_data = get_completed_step_data(state_dir, sample, assembler, config, step)
        if step_data and "path" in step_data and "result" in step_data:
            json_file = step_data["path"]
            if os.path.exists(json_file):
                logger.info(f"✅ Resuming: BUSCO {lineage} already completed for {sample}/{assembler}/{config}")
                return json_file
    
    os.makedirs(out_dir, exist_ok=True)
    out_name = f"{sample}_{assembler}_{config}_{lineage}"
    
    # Check if the BUSCO output directory already exists and has completed JSON
    busco_subdir = os.path.join(out_dir, out_name)
    json_candidates = glob.glob(os.path.join(busco_subdir, "short_summary.*.json"))
    
    if json_candidates and not force:
        logger.info(f"✅ Resuming: Found existing BUSCO {lineage} results for {sample}/{assembler}/{config}")
        json_file = json_candidates[0]
    else:
        # Run BUSCO
        cmd = (
            f"busco -i {fasta_file} -o {out_name} -l {lineage} -m genome "
            f"--cpu {threads} --out_path {out_dir} --force"
        )
        ok, _ = run_cmd(cmd, desc=f"BUSCO {lineage} for {sample}/{assembler}/{config}")
        if not ok:
            if state_dir:
                update_state(state_dir, sample, assembler, config, step, "failed", 
                             error=f"busco_{lineage}_failed")
            return None

        # Check for JSON after running
        json_candidates = glob.glob(os.path.join(busco_subdir, "short_summary.*.json"))
        if not json_candidates:
            logger.error(f"No BUSCO JSON found for lineage={lineage}, {sample}/{assembler}/{config}")
            if state_dir:
                update_state(state_dir, sample, assembler, config, step, "failed", 
                             error=f"busco_{lineage}_missing_json")
            return None
        json_file = json_candidates[0]
    
    # Parse the BUSCO results
    busco_results = parse_busco_json(json_file)
    
    # Update state
    if state_dir:
        update_state(
            state_dir, sample, assembler, config, step, "completed", 
            path=json_file, 
            result=busco_results,
            timestamp=datetime.now().isoformat()
        )
    
    return json_file

# --------------------------------------------------------------------
# 7. Final Summaries
# --------------------------------------------------------------------
def create_comparison_summary(results, out_dir):
    """
    Make a DataFrame from 'results' (list of dicts), then sort by:
      busco_mollusca desc, busco_eukaryota desc, n50 desc
    Save to CSV and a markdown summary.
    """
    df = pd.DataFrame(results)
    if df.empty:
        logger.warning("No results to summarize.")
        return

    # Sort
    df = df.sort_values(
        by=["busco_mollusca", "busco_eukaryota", "n50"], 
        ascending=False
    )

    best = df.iloc[0]

    # Save CSV
    csv_path = os.path.join(out_dir, "assembly_comparison_summary.csv")
    df.to_csv(csv_path, index=False)

    # Save Markdown
    md_path = os.path.join(out_dir, "assembly_comparison_summary.md")
    with open(md_path, 'w') as f:
        f.write("# Long-read Assembly Comparison Summary\n\n")
        f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write("Sorted by **Mollusca completeness** desc, then **Eukaryota completeness** desc, then **N50** desc.\n\n")

        f.write(df.to_markdown(index=False))
        f.write("\n\n")
        f.write("## Best Assembly\n\n")
        f.write(f"- Sample: **{best['sample']}**\n")
        f.write(f"- Assembler: **{best['assembler']}**\n")
        f.write(f"- Config: **{best['config']}**\n")
        f.write(f"- BUSCO (Mollusca) = {best['busco_mollusca']:.2f}%\n")
        f.write(f"- BUSCO (Eukaryota) = {best['busco_eukaryota']:.2f}%\n")
        f.write(f"- N50 = {best['n50']:,}\n")
        f.write(f"- GC% = {best['gc']:.2f}\n")

    logger.info(f"Comparison summary saved:\n  {csv_path}\n  {md_path}")

def collect_results_from_state(output_dir):
    """
    Collect results from the state file to recover previous runs.
    Returns a list of result dictionaries.
    """
    state = read_state(output_dir)
    results = []
    
    for sample, sample_data in state.items():
        if "assemblers" not in sample_data:
            continue
        
        for asm_key, asm_data in sample_data["assemblers"].items():
            # Parse assembler and config from asm_key (e.g., "FLYE_pacbio-hifi")
            parts = asm_key.split("_", 1)
            if len(parts) < 2:
                continue
            assembler, config = parts
            
            # Check if assembly, busco, and gc steps are completed
            assembly_data = asm_data.get("assembly", {})
            busco_mollusca_data = asm_data.get("busco_mollusca", {})
            busco_eukaryota_data = asm_data.get("busco_eukaryota", {})
            gc_data = asm_data.get("gc", {})
            
            if assembly_data.get("status") != "completed":
                continue
            
            # Gather results
            row = {
                "sample": sample,
                "assembler": assembler,
                "config": config,
            }
            
            # Add BUSCO mollusca results if available
            if busco_mollusca_data.get("status") == "completed" and "result" in busco_mollusca_data:
                result = busco_mollusca_data["result"]
                row["busco_mollusca"] = result.get("complete_percent", 0.0)
                row["missing_mollusca"] = result.get("missing_percent", 0.0)
                row["n50"] = result.get("n50", 0)
                row["num_contigs"] = result.get("num_contigs", 0)
                row["total_length"] = result.get("total_length", 0)
            else:
                row["busco_mollusca"] = 0.0
                row["missing_mollusca"] = 0.0
                row["n50"] = 0
                row["num_contigs"] = 0
                row["total_length"] = 0
            
            # Add BUSCO eukaryota results if available
            if busco_eukaryota_data.get("status") == "completed" and "result" in busco_eukaryota_data:
                result = busco_eukaryota_data["result"]
                row["busco_eukaryota"] = result.get("complete_percent", 0.0)
                row["missing_eukaryota"] = result.get("missing_percent", 0.0)
                # Use eukaryota values if mollusca is missing or has lower values
                if row["n50"] == 0:
                    row["n50"] = result.get("n50", 0)
                if row["num_contigs"] == 0:
                    row["num_contigs"] = result.get("num_contigs", 0)
                if row["total_length"] == 0:
                    row["total_length"] = result.get("total_length", 0)
            else:
                row["busco_eukaryota"] = 0.0
                row["missing_eukaryota"] = 0.0
            
            # Add GC content if available
            if gc_data.get("status") == "completed" and "value" in gc_data:
                row["gc"] = gc_data["value"]
            else:
                row["gc"] = 0.0
            
            results.append(row)
    
    return results

def check_dependencies():
    """
    Check for required command-line tools and Python packages.
    Exits if any are missing.
    """
    import shutil
    import importlib

    required_tools = [
        "flye",
        "canu",
        "busco",
        "seqkit"
    ]

    required_python_packages = [
        "pandas"
    ]

    missing_tools = [tool for tool in required_tools if shutil.which(tool) is None]
    missing_packages = [pkg for pkg in required_python_packages if importlib.util.find_spec(pkg) is None]

    if missing_tools or missing_packages:
        print("\n❌ Missing dependencies:")
        if missing_tools:
            print("🔧 Missing tools:")
            for tool in missing_tools:
                print(f"  - {tool}")
        if missing_packages:
            print("🐍 Missing Python packages:")
            for pkg in missing_packages:
                print(f"  - {pkg}")
        print("\n💡 Please install the missing components and try again.")
        print("💡 If you're using Conda/Mamba, you can run:")
        print("    mamba install -c bioconda flye canu busco seqkit")
        print()
        sys.exit(1)

    print("✅ All dependencies satisfied.\n")

# --------------------------------------------------------------------
# 8. Main
# --------------------------------------------------------------------
def main():
    check_dependencies()
    parser = argparse.ArgumentParser(
        description="Long-read assembly pipeline comparing FLYE and CANU with BUSCO and GC%. With resume functionality."
    )
    parser.add_argument("-i", "--input_dir", required=True, help="Directory with long read files (FASTQ/FASTA)")
    parser.add_argument("-o", "--output_dir", default="Longread_Assembly_Comparison", help="Output directory")
    parser.add_argument("-t", "--threads", type=int, default=16, help="Number of CPU threads")
    parser.add_argument("-g", "--genome_size", required=True, help="Approximate genome size (e.g., '50m' for 50 Mbp)")
    parser.add_argument("--tech", default="pacbio-hifi", choices=["pacbio-hifi", "pacbio-raw", "nano-raw", "nano-hq", "nano-corr"], 
                        help="Long read technology type")
    parser.add_argument("--clean-after", action="store_true", help="Delete intermediate build folders after completion")
    parser.add_argument("--force", action="store_true", help="Force re-run all steps, ignoring previous results")
    parser.add_argument("--force-assembly", action="store_true", help="Force re-run assembly steps only")
    parser.add_argument("--force-busco", action="store_true", help="Force re-run BUSCO steps only")
    parser.add_argument("--force-gc", action="store_true", help="Force re-run GC calculation only")
    parser.add_argument("--show-state", action="store_true", help="Show current state information and exit")
    parser.add_argument("--reset-state", action="store_true", help="Reset pipeline state before running")
    args = parser.parse_args()

    log_fp = setup_logging(args.output_dir)
    logger.info(f"=== Long-read Assembly Pipeline: FLYE + CANU ===")
    logger.info(f"Input dir: {args.input_dir}")
    logger.info(f"Output dir: {args.output_dir}")
    logger.info(f"Technology: {args.tech}")
    logger.info(f"Genome size: {args.genome_size}")
    logger.info(f"Threads: {args.threads}")
    logger.info(f"Logs -> {log_fp}")
    
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
            for sample, data in state.items():
                logger.info(f"Sample: {sample}")
                for asm_key, asm_data in data.get("assemblers", {}).items():
                    logger.info(f"  {asm_key}:")
                    for step, step_data in asm_data.items():
                        status = step_data.get("status", "unknown")
                        status_symbol = "✅" if status == "completed" else "❌"
                        logger.info(f"    {step}: {status_symbol} {status}")
            
            # Print summary of completed steps
            completed_assemblies = sum(
                1 for s in state.values() 
                for a in s.get("assemblers", {}).values() 
                if a.get("assembly", {}).get("status") == "completed"
            )
            completed_busco = sum(
                1 for s in state.values() 
                for a in s.get("assemblers", {}).values() 
                for step in ["busco_mollusca", "busco_eukaryota"]
                if a.get(step, {}).get("status") == "completed"
            )
            completed_gc = sum(
                1 for s in state.values() 
                for a in s.get("assemblers", {}).values() 
                if a.get("gc", {}).get("status") == "completed"
            )
            
            logger.info(f"\nSummary: {completed_assemblies} assemblies, {completed_busco} BUSCO runs, {completed_gc} GC calculations completed")
        
        # Collect and show result summary if available
        results = collect_results_from_state(args.output_dir)
        if results:
            df = pd.DataFrame(results)
            df = df.sort_values(by=["busco_mollusca", "busco_eukaryota", "n50"], ascending=False)
            logger.info("\nResults summary (top 5):")
            pd.set_option('display.max_columns', None)
            logger.info(df.head(5).to_string())
        
        sys.exit(0)
    
    # Create subdirs
    asm_dir   = os.path.join(args.output_dir, "assemblies")
    busco_dir = os.path.join(args.output_dir, "busco")
    os.makedirs(asm_dir, exist_ok=True)
    os.makedirs(busco_dir, exist_ok=True)

    # Find long read files
    long_read_files = []
    for ext in [".fastq", ".fq", ".fasta", ".fa", ".fastq.gz", ".fq.gz", ".fasta.gz", ".fa.gz"]:
        long_read_files.extend(glob.glob(os.path.join(args.input_dir, f"*{ext}")))
    
    if not long_read_files:
        logger.error("No long read files found. Exiting.")
        sys.exit(1)

    # Map FLYE read type to CANU read type
    flye_to_canu_map = {
        "pacbio-hifi": "pacbio-hifi",
        "pacbio-raw": "pacbio-raw",
        "nano-raw": "nanopore-raw",
        "nano-hq": "nanopore-raw",  # Use raw for HQ nanopore reads
        "nano-corr": "nanopore-corr"
    }
    
    # Determine which steps to force based on command-line args
    force_assembly = args.force or args.force_assembly
    force_busco = args.force or args.force_busco
    force_gc = args.force or args.force_gc
    
    # Detect if we are resuming by checking the state file
    resuming = os.path.exists(state_file_path(args.output_dir)) and not args.force and not args.reset_state
    if resuming:
        logger.info("🔄 Resuming pipeline from previous run")
        # Load the existing results from state
        existing_results = collect_results_from_state(args.output_dir)
        # Create a set of already completed sample+assembler+config combinations
        completed_configs = set(
            f"{r['sample']}_{r['assembler']}_{r['config']}" 
            for r in existing_results
        )
        logger.info(f"Found {len(existing_results)} completed assemblies from previous runs")
    else:
        completed_configs = set()
        existing_results = []
    
    # We'll store final results in a list of dicts
    all_results = list(existing_results)
    lineages = ["mollusca_odb10", "eukaryota_odb10"]  # run both lineages

    for reads_file in long_read_files:
        sample_name = os.path.basename(reads_file).split(".")[0]
        logger.info(f"\n=== Sample: {sample_name} ===")

        # FLYE assembly
        flye_config = args.tech
        config_key = f"{sample_name}_FLYE_{flye_config}"
        
        if config_key not in completed_configs or force_assembly:
            flye_asm = run_flye_assembly(
                reads=reads_file,
                output_dir=asm_dir,
                sample_name=sample_name,
                read_type=flye_config,
                genome_size=args.genome_size,
                threads=args.threads,
                state_dir=args.output_dir,
                force=force_assembly
            )
            
            if flye_asm:
                # Get GC content
                gc_flye = get_gc_content(
                    flye_asm, 
                    output_dir=args.output_dir, 
                    sample=sample_name, 
                    assembler="FLYE", 
                    config=flye_config, 
                    force=force_gc
                )
                
                # Run BUSCO for both lineages
                busco_data = {}
                for lin in lineages:
                    busco_out_dir = os.path.join(busco_dir, f"flye_{flye_config}_{lin}")
                    busco_json_file = run_busco(
                        flye_asm, busco_out_dir, lin, args.threads, 
                        sample_name, "FLYE", flye_config,
                        state_dir=args.output_dir, force=force_busco
                    )
                    if busco_json_file:
                        busco_data[lin] = parse_busco_json(busco_json_file)
                
                # Add row to results if not already present or if forced
                if config_key not in completed_configs or force_assembly or force_busco or force_gc:
                    moll = busco_data.get("mollusca_odb10", {})
                    euk  = busco_data.get("eukaryota_odb10", {})
                    row = {
                        "sample": sample_name,
                        "assembler": "FLYE",
                        "config": flye_config,
                        "busco_mollusca": moll.get("complete_percent", 0.0),
                        "missing_mollusca": moll.get("missing_percent", 0.0),
                        "busco_eukaryota": euk.get("complete_percent", 0.0),
                        "missing_eukaryota": euk.get("missing_percent", 0.0),
                        "n50": max(moll.get("n50",0), euk.get("n50",0)),
                        "num_contigs": max(moll.get("num_contigs",0), euk.get("num_contigs",0)),
                        "total_length": max(moll.get("total_length",0), euk.get("total_length",0)),
                        "gc": gc_flye
                    }
                    # Remove existing entry if present (for updates)
                    all_results = [r for r in all_results if not (r["sample"] == sample_name and 
                                                                r["assembler"] == "FLYE" and 
                                                                r["config"] == flye_config)]
                    all_results.append(row)
        else:
            logger.info(f"✅ Skipping FLYE {flye_config} for {sample_name} (already completed)")
        
        # CANU assembly
        canu_config = flye_to_canu_map.get(args.tech, "pacbio-hifi")
        config_key = f"{sample_name}_CANU_{canu_config}"
        
        if config_key not in completed_configs or force_assembly:
            canu_asm = run_canu_assembly(
                reads=reads_file,
                output_dir=asm_dir,
                sample_name=sample_name,
                read_type=canu_config,
                genome_size=args.genome_size,
                threads=args.threads,
                state_dir=args.output_dir,
                force=force_assembly
            )
            
            if canu_asm:
                # Get GC content
                gc_canu = get_gc_content(
                    canu_asm, 
                    output_dir=args.output_dir, 
                    sample=sample_name, 
                    assembler="CANU", 
                    config=canu_config, 
                    force=force_gc
                )
                
                # Run BUSCO for both lineages
                busco_data = {}
                for lin in lineages:
                    busco_out_dir = os.path.join(busco_dir, f"canu_{canu_config}_{lin}")
                    busco_json_file = run_busco(
                        canu_asm, busco_out_dir, lin, args.threads, 
                        sample_name, "CANU", canu_config,
                        state_dir=args.output_dir, force=force_busco
                    )
                    if busco_json_file:
                        busco_data[lin] = parse_busco_json(busco_json_file)
                
                # Add row to results
                if config_key not in completed_configs or force_assembly or force_busco or force_gc:
                    moll = busco_data.get("mollusca_odb10", {})
                    euk  = busco_data.get("eukaryota_odb10", {})
                    row = {
                        "sample": sample_name,
                        "assembler": "CANU",
                        "config": canu_config,
                        "busco_mollusca": moll.get("complete_percent", 0.0),
                        "missing_mollusca": moll.get("missing_percent", 0.0),
                        "busco_eukaryota": euk.get("complete_percent", 0.0),
                        "missing_eukaryota": euk.get("missing_percent", 0.0),
                        "n50": max(moll.get("n50",0), euk.get("n50",0)),
                        "num_contigs": max(moll.get("num_contigs",0), euk.get("num_contigs",0)),
                        "total_length": max(moll.get("total_length",0), euk.get("total_length",0)),
                        "gc": gc_canu
                    }
                    # Remove existing entry if present (for updates)
                    all_results = [r for r in all_results if not (r["sample"] == sample_name and 
                                                                r["assembler"] == "CANU" and 
                                                                r["config"] == canu_config)]
                    all_results.append(row)
        else:
            logger.info(f"✅ Skipping CANU {canu_config} for {sample_name} (already completed)")

    # Create final summary
    create_comparison_summary(all_results, args.output_dir)
    
    # Clean up if requested
    if args.clean_after:
        for reads_file in long_read_files:
            sample_name = os.path.basename(reads_file).split(".")[0]
            cleanup_dirs = [
                os.path.join(asm_dir, f"{sample_name}_flye_{args.tech}"),
                os.path.join(asm_dir, f"{sample_name}_canu_{flye_to_canu_map.get(args.tech)}")
            ]
            for d in cleanup_dirs:
                if os.path.exists(d):
                    shutil.rmtree(d)
                    logger.info(f"🧹 Cleaned up temp folder: {d}")
    
    logger.info("Pipeline complete.")

if __name__ == "__main__":
    main()