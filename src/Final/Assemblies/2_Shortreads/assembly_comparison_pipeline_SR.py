#!/usr/bin/env python3
"""
All-vs-All Assembly Pipeline (SOAP + GapCloser, MaSuRCA permissive/strict, Captus)
with BUSCO (mollusca_odb10 & eukaryota_odb10) and GC% via seqkit.
With resume functionality to restart from failed steps.

Author: Max Borgmann
Date: 03.04.2025

This version uses separate conda environments (via mamba) to avoid conflicts:
  - soap_env
  - masurca_env
  - captus_env

You only need to run this script once; it will:
  1) Create the above 3 envs if missing.
  2) Run the entire pipeline in each env as needed.
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
from datetime import datetime
import hashlib
import tempfile
import time

# --------------------------------------------------------------------------------
# 0. Create/Check Conda Environments via Mamba
# --------------------------------------------------------------------------------

def get_existing_envs_via_mamba():
    """
    Return a set of all environment *names* recognized by 'mamba env list --json'.
    """
    cmd = "mamba env list --json"
    existing = set()
    try:
        result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
        data = json.loads(result.stdout)
        for path in data["envs"]:
            # The env name is typically the basename of the path
            env_name = os.path.basename(path)
            existing.add(env_name)
    except Exception as e:
        print(f"Warning: could not parse existing envs from 'mamba env list': {e}")
    return existing

def setup_conda_envs():
    """
    Creates three envs via mamba if they do not already exist:
      - soap_env
      - masurca_env
      - captus_env
    Each is pinned to python=3.10 and includes busco=5.8.2 & seqkit=2.3.1.
    """
    existing_envs = get_existing_envs_via_mamba()

    # --- soap_env ---
    if "soap_env" not in existing_envs:
        print("Creating conda environment: soap_env")
        cmd = (
            "mamba create -y -n soap_env "
            "-c conda-forge -c bioconda "
            "python=3.10 soapdenovo2-gapcloser=2.04 "
            "busco=5.8.2 seqkit=2.3.1"
        )
        subprocess.run(cmd, shell=True, check=True)
    else:
        print("Environment already exists: soap_env")

    # --- masurca_env ---
    if "masurca_env" not in existing_envs:
        print("Creating conda environment: masurca_env")
        cmd = (
            "mamba create -y -n masurca_env "
            "-c conda-forge -c bioconda "
            "python=3.10 masurca=4.0.9 "
            "busco=5.8.2 seqkit=2.3.1"
        )
        subprocess.run(cmd, shell=True, check=True)
    else:
        print("Environment already exists: masurca_env")

    # --- captus_env ---
    if "captus_env" not in existing_envs:
        print("Creating conda environment: captus_env")
        cmd = (
            "mamba create -y -n captus_env "
            "-c conda-forge -c bioconda "
            "python=3.10 captus=1.3.0 "
            "busco=5.8.2 seqkit=2.3.1"
        )
        subprocess.run(cmd, shell=True, check=True)
    else:
        print("Environment already exists: captus_env")

def run_in_env(env_name, cmd, cwd=None, desc=None):
    """
    Run 'cmd' inside conda environment 'env_name' with mamba run (or conda run).
    Returns (ok, stdout).
    """
    if desc:
        logging.info(f"[{env_name}] {desc}")
    full_cmd = f"mamba run -n {env_name} {cmd}"
    logging.info(f"Running in {env_name}:\n  {full_cmd}")
    try:
        result = subprocess.run(full_cmd, shell=True, check=True, cwd=cwd,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if result.stdout.strip():
            logging.debug(f"STDOUT: {result.stdout.strip()}")
        if result.stderr.strip():
            logging.debug(f"STDERR: {result.stderr.strip()}")
        return (True, result.stdout)
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed in env={env_name}: {cmd}")
        logging.error(f"Return code: {e.returncode}")
        logging.error(f"STDOUT: {e.stdout}")
        logging.error(f"STDERR: {e.stderr}")
        return (False, e.stderr)

# --------------------------------------------------------------------------------
# 1. Logging Setup
# --------------------------------------------------------------------------------
logger = logging.getLogger("assembly_comparison")
logger.setLevel(logging.DEBUG)

def setup_logging(log_dir, log_name="assembly_comparison.log"):
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

# --------------------------------------------------------------------------------
# 2. State Management for Resume Functionality
# --------------------------------------------------------------------------------
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
    """Read the current pipeline state from JSON, or {} if not found/corrupted."""
    sf = state_file_path(output_dir)
    if not os.path.exists(sf):
        return {}
    try:
        with open(sf, 'r') as f:
            return json.load(f)
    except (json.JSONDecodeError, FileNotFoundError):
        logger.warning("State file corrupted or not found. Starting fresh state.")
        return {}

def write_state(state, output_dir):
    """Write the pipeline state to disk with atomic operations."""
    sf = state_file_path(output_dir)
    os.makedirs(output_dir, exist_ok=True)
    with tempfile.NamedTemporaryFile('w', dir=output_dir, delete=False) as tmp:
        json.dump(state, tmp, indent=2)
        tmp_path = tmp.name
    shutil.move(tmp_path, sf)
    logger.debug(f"Updated pipeline state file: {sf}")

def update_state(output_dir, sample, assembler, config, step, status, **kwargs):
    """
    Update the pipeline state for a specific step.
    step can be: "assembly", "busco_mollusca", "busco_eukaryota", "gc", etc.
    """
    state = read_state(output_dir)
    if sample not in state:
        state[sample] = {"last_update": datetime.now().isoformat(), "assemblers": {}}
    asm_key = f"{assembler}_{config}"
    if asm_key not in state[sample]["assemblers"]:
        state[sample]["assemblers"][asm_key] = {}

    state[sample]["assemblers"][asm_key][step] = {"status": status, **kwargs}
    state[sample]["last_update"] = datetime.now().isoformat()
    write_state(state, output_dir)

def is_step_completed(output_dir, sample, assembler, config, step, verify_file=None):
    """
    Check if a step has been completed successfully.
    If verify_file is provided, also verify that file's existence/hash.
    """
    st = read_state(output_dir)
    try:
        step_state = st[sample]["assemblers"][f"{assembler}_{config}"][step]
        if step_state["status"] != "completed":
            return False
        if verify_file and "path" in step_state:
            stored_path = step_state["path"]
            if not os.path.exists(stored_path):
                logger.warning(f"File marked as completed but not found: {stored_path}")
                return False
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
    st = read_state(output_dir)
    try:
        return st[sample]["assemblers"][f"{assembler}_{config}"][step]
    except (KeyError, TypeError):
        return None

# --------------------------------------------------------------------------------
# 3. Utility: parse BUSCO JSON, get seqkit stats (run_in_env)
# --------------------------------------------------------------------------------
def parse_busco_json(json_file):
    """
    Parse BUSCO v5+ JSON to extract:
      - 'complete_percent'
      - 'missing_percent'
      - 'n50'
      - 'num_contigs'
      - 'total_length'
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

def get_gc_content(env_name, fasta_file, output_dir=None, sample=None, assembler=None, config=None, force=False):
    """
    Use seqkit to get GC% of the entire assembly, with resume.
    We run seqkit inside the given environment.
    """
    if not fasta_file or not os.path.exists(fasta_file):
        return 0.0

    # Check state
    if output_dir and sample and assembler and config and not force:
        step_data = get_completed_step_data(output_dir, sample, assembler, config, "gc")
        if step_data and "value" in step_data:
            logger.info(f"✅ Resuming: GC content already known for {sample}/{assembler}/{config}")
            return step_data["value"]

    out_file = f"{fasta_file}.seqkit_stats.txt"
    cmd = f"seqkit stats -a {fasta_file} -o {out_file}"
    ok, _ = run_in_env(env_name, cmd, desc="seqkit stats")
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

# --------------------------------------------------------------------------------
# 4. SOAP with multi-step + GapCloser
# --------------------------------------------------------------------------------
def run_soapdenovo2_with_gapcloser(r1, r2, output_dir, sample_name, k_value, threads=16, state_dir=None, force=False):
    """
    Run SOAPdenovo2 (with the soap_env). Then run GapCloser (also in soap_env).
    """
    env_name = "soap_env"
    assembler = "SOAP"
    config = f"k{k_value}"

    # Check if assembly done
    if not force and state_dir and is_step_completed(state_dir, sample_name, assembler, config, "assembly"):
        step_data = get_completed_step_data(state_dir, sample_name, assembler, config, "assembly")
        if step_data and "path" in step_data:
            asm_path = step_data["path"]
            if os.path.exists(asm_path) and os.path.getsize(asm_path) > 0:
                logger.info(f"✅ Resuming: SOAP assembly k={k_value} for {sample_name}")
                return asm_path

    os.makedirs(output_dir, exist_ok=True)
    soap_dir = os.path.join(output_dir, f"{sample_name}_soap_k{k_value}")
    os.makedirs(soap_dir, exist_ok=True)

    config_file = os.path.join(soap_dir, "soap_config.txt")
    with open(config_file, 'w') as f:
        f.write(f"""max_rd_len=150
[LIB]
avg_ins=150
reverse_seq=0
asm_flags=3
rank=1
pair_num_cutoff=5
map_len=40
q1={r1}
q2={r2}
""")

    prefix = os.path.join(soap_dir, f"{sample_name}_soap_k{k_value}")
    pregraph_file = f"{prefix}.preGraphBasic"
    contig_file   = f"{prefix}.contig"
    map_file      = f"{prefix}.readOnContig.readInGap.peGrads"
    scaff_file    = f"{prefix}.scafSeq"
    gapclosed_file= f"{prefix}.gapclosed.fasta"

    # 1) pregraph
    if not os.path.exists(pregraph_file) or force:
        pre_cmd = f"SOAPdenovo-{k_value}mer pregraph -s {config_file} -K {k_value} -p {threads} -o {prefix}"
        ok, _ = run_in_env(env_name, pre_cmd, desc=f"SOAP pregraph k={k_value}")
        if not ok:
            if state_dir:
                update_state(state_dir, sample_name, assembler, config, "assembly", "failed", error="pregraph_failed")
            return None
    else:
        logger.info(f"✅ Skipping SOAP pregraph k={k_value}; already done.")

    # 2) contig
    if not os.path.exists(contig_file) or force:
        contig_cmd = f"SOAPdenovo-{k_value}mer contig -g {prefix}"
        ok, _ = run_in_env(env_name, contig_cmd, desc=f"SOAP contig k={k_value}")
        if not ok:
            if state_dir:
                update_state(state_dir, sample_name, assembler, config, "assembly", "failed", error="contig_failed")
            return None
    else:
        logger.info(f"✅ Skipping SOAP contig k={k_value}; already done.")

    # 3) map
    if not os.path.exists(map_file) or force:
        map_cmd = f"SOAPdenovo-{k_value}mer map -s {config_file} -g {prefix} -p {threads}"
        ok, _ = run_in_env(env_name, map_cmd, desc=f"SOAP map k={k_value}")
        if not ok:
            if state_dir:
                update_state(state_dir, sample_name, assembler, config, "assembly", "failed", error="map_failed")
            return None
    else:
        logger.info(f"✅ Skipping SOAP map k={k_value}; already done.")

    # 4) scaff
    if not os.path.exists(scaff_file) or force:
        scaff_cmd = f"SOAPdenovo-{k_value}mer scaff -g {prefix} -F"
        ok, _ = run_in_env(env_name, scaff_cmd, desc=f"SOAP scaff k={k_value}")
        if not ok:
            if state_dir:
                update_state(state_dir, sample_name, assembler, config, "assembly", "failed", error="scaff_failed")
            return None
    else:
        logger.info(f"✅ Skipping SOAP scaff k={k_value}; already done.")

    # 5) GapCloser
    if not os.path.exists(gapclosed_file) or force or os.path.getsize(gapclosed_file) == 0:
        gapcloser_cmd = f"GapCloser -b {config_file} -a {prefix}.scafSeq -o {gapclosed_file} -t {threads}"
        ok, _ = run_in_env(env_name, gapcloser_cmd, desc="GapCloser")
        if not ok:
            if state_dir:
                update_state(state_dir, sample_name, assembler, config, "assembly", "failed", error="gapcloser_failed")
            return None
    else:
        logger.info("✅ Skipping GapCloser; already done or file present.")

    if not os.path.exists(gapclosed_file):
        logger.error(f"Gapclosed file not found for SOAP k={k_value}")
        if state_dir:
            update_state(state_dir, sample_name, assembler, config, "assembly", "failed", error="gapclosed_missing")
        return None

    final_asm = os.path.join(output_dir, f"{sample_name}_soap_k{k_value}.fasta")
    shutil.copy2(gapclosed_file, final_asm)

    # Update state
    if state_dir:
        update_state(
            state_dir, sample_name, assembler, config, "assembly", "completed",
            path=final_asm,
            hash=get_file_hash(final_asm),
            timestamp=datetime.now().isoformat()
        )
    return final_asm

# --------------------------------------------------------------------------------
# 5. MaSuRCA (Permissive vs Stringent)
# --------------------------------------------------------------------------------
def run_masurca(r1, r2, output_dir, sample_name, config_name, threads=16, state_dir=None, force=False):
    """
    Run MaSuRCA in masurca_env with a custom config.
    'permissive' vs 'stringent'
    """
    env_name = "masurca_env"
    assembler = "MaSuRCA"
    config = config_name

    if not force and state_dir and is_step_completed(state_dir, sample_name, assembler, config, "assembly"):
        step_data = get_completed_step_data(state_dir, sample_name, assembler, config, "assembly")
        if step_data and "path" in step_data:
            asm_path = step_data["path"]
            if os.path.exists(asm_path) and os.path.getsize(asm_path) > 0:
                logger.info(f"✅ Resuming: MaSuRCA {config_name} for {sample_name}")
                return asm_path

    os.makedirs(output_dir, exist_ok=True)
    ms_dir = os.path.join(output_dir, f"{sample_name}_masurca_{config_name}")
    os.makedirs(ms_dir, exist_ok=True)

    if config_name == "permissive":
        kmer_size = 31
        kmer_count= 2
        cgw_error = 0.15
        ovl_mer   = 30
        do_trim   = 0
    else:
        kmer_size = 61
        kmer_count= 3
        cgw_error = 0.10
        ovl_mer   = 45
        do_trim   = 1

    jf_size = 1300000000
    config_file = os.path.join(ms_dir, "masurca_config.txt")
    with open(config_file, 'w') as f:
        f.write(f"""# MaSuRCA configuration
DATA
PE = pe 150 15 {r1} {r2}
END

PARAMETERS
GRAPH_KMER_SIZE = {kmer_size}
USE_LINKING_MATES = 1
LIMIT_JUMP_COVERAGE = 300
CA_PARAMETERS = ovlMerSize={ovl_mer} cgwErrorRate={cgw_error} ovlHashBits=25 ovlHashBlockLength=100000000 utgMemory=650GB obtMemory=650GB
KMER_COUNT_THRESHOLD = {kmer_count}
NUM_THREADS = {threads}
JF_SIZE = {jf_size}
DO_HOMOPOLYMER_TRIM = {do_trim}
CLOSE_GAPS = 1
SOAP_ASSEMBLY = 0
END
""")

    final_candidates = [
        os.path.join(ms_dir, "CA", "final.genome.scf.fasta"),
        os.path.join(ms_dir, "CA", "9-terminator", "genome.scf.fasta")
    ]
    found_existing = False
    final_asm_path = None

    for cand in final_candidates:
        if os.path.exists(cand) and os.path.getsize(cand) > 0 and not force:
            found_existing = True
            final_asm_path = cand
            logger.info(f"✅ Resuming: Found existing MaSuRCA assembly: {cand}")
            break

    if not found_existing or force:
        # 1) generate assemble.sh
        assemble_sh = os.path.join(ms_dir, "assemble.sh")
        if not os.path.exists(assemble_sh) or force:
            cmd1 = f"masurca {config_file}"
            ok, _ = run_in_env(env_name, cmd1, cwd=ms_dir, desc="MaSuRCA config->assemble.sh")
            if not ok:
                if state_dir:
                    update_state(state_dir, sample_name, assembler, config, "assembly", "failed", error="masurca_config_failed")
                return None
        else:
            logger.info(f"✅ Skipping MaSuRCA config step; assemble.sh already exists.")

        # 2) run assemble.sh
        cmd2 = "./assemble.sh"
        ok, _ = run_in_env(env_name, cmd2, cwd=ms_dir, desc=f"MaSuRCA run {config_name}")
        if not ok:
            if state_dir:
                update_state(state_dir, sample_name, assembler, config, "assembly", "failed", error="masurca_assembly_failed")
            return None

        # check final scf
        for cand in final_candidates:
            if os.path.exists(cand) and os.path.getsize(cand) > 0:
                final_asm_path = cand
                break

        if not final_asm_path:
            logger.error(f"No final scf found for {config_name} in {ms_dir}")
            if state_dir:
                update_state(state_dir, sample_name, assembler, config, "assembly", "failed", error="masurca_missing_final")
            return None

    out_asm = os.path.join(output_dir, f"{sample_name}_masurca_{config_name}.fasta")
    shutil.copy2(final_asm_path, out_asm)

    # Update state
    if state_dir:
        update_state(
            state_dir, sample_name, assembler, config, "assembly", "completed",
            path=out_asm,
            hash=get_file_hash(out_asm),
            timestamp=datetime.now().isoformat()
        )
    return out_asm

# --------------------------------------------------------------------------------
# 6. Captus
# --------------------------------------------------------------------------------
def run_captus_assembly(r1, r2, output_dir, sample_name, threads=16, state_dir=None, force=False):
    """
    Run Captus with WGS preset in 'captus_env'.
    """
    env_name = "captus_env"
    assembler = "Captus"
    config = "WGS"

    if not force and state_dir and is_step_completed(state_dir, sample_name, assembler, config, "assembly"):
        step_data = get_completed_step_data(state_dir, sample_name, assembler, config, "assembly")
        if step_data and "path" in step_data:
            asm_path = step_data["path"]
            if os.path.exists(asm_path) and os.path.getsize(asm_path) > 0:
                logger.info(f"✅ Resuming: Captus WGS for {sample_name}")
                return asm_path

    os.makedirs(output_dir, exist_ok=True)
    captus_dir = os.path.join(output_dir, f"{sample_name}_captus_wgs_run")
    os.makedirs(captus_dir, exist_ok=True)

    final_asm_path = os.path.join(captus_dir, f"{sample_name}__captus-asm", "assembly.fasta")
    if os.path.exists(final_asm_path) and os.path.getsize(final_asm_path) > 0 and not force:
        logger.info(f"✅ Resuming: Found existing Captus assembly: {final_asm_path}")
    else:
        clean_dir = os.path.join(captus_dir, "01_clean_reads")
        os.makedirs(clean_dir, exist_ok=True)

        r1_link = os.path.join(clean_dir, f"{sample_name}_R1.fastq")
        r2_link = os.path.join(clean_dir, f"{sample_name}_R2.fastq")
        for ln in [r1_link, r2_link]:
            if os.path.islink(ln):
                os.unlink(ln)
        os.symlink(os.path.abspath(r1), r1_link)
        os.symlink(os.path.abspath(r2), r2_link)

        cmd = (
            f"captus assemble -r {clean_dir} -o {captus_dir} "
            f"--threads {threads} --ram 32 --preset WGS "
            f"--k-list 31,39,49,69,89,109,129,149,169 --min-count 3 "
            f"--prune-level 2 --no-mercy"
        )
        ok, _ = run_in_env(env_name, cmd, desc="Captus assembly")
        if not ok:
            if state_dir:
                update_state(state_dir, sample_name, assembler, config, "assembly", "failed", error="captus_failed")
            return None

        if not os.path.exists(final_asm_path):
            logger.error("Captus final assembly not found.")
            if state_dir:
                update_state(state_dir, sample_name, assembler, config, "assembly", "failed", error="captus_missing_final")
            return None

    out_fa = os.path.join(output_dir, f"{sample_name}_captus_wgs.fasta")
    shutil.copy2(final_asm_path, out_fa)

    # Update state
    if state_dir:
        update_state(
            state_dir, sample_name, assembler, config, "assembly", "completed",
            path=out_fa,
            hash=get_file_hash(out_fa),
            timestamp=datetime.now().isoformat()
        )
    return out_fa

# --------------------------------------------------------------------------------
# 7. BUSCO
# --------------------------------------------------------------------------------
def run_busco(env_name, fasta_file, out_dir, lineage, threads, sample, assembler, config, state_dir=None, force=False):
    """
    Runs BUSCO in env=env_name for 'fasta_file' with the given lineage, storing JSON in out_dir.
    Updates pipeline state on success.
    """
    step = f"busco_{lineage.split('_')[0]}"  # e.g. "busco_mollusca"

    if not force and state_dir and is_step_completed(state_dir, sample, assembler, config, step):
        step_data = get_completed_step_data(state_dir, sample, assembler, config, step)
        if step_data and "path" in step_data and "result" in step_data:
            json_file = step_data["path"]
            if os.path.exists(json_file):
                logger.info(f"✅ Resuming: BUSCO {lineage} for {sample}/{assembler}/{config}")
                return json_file

    os.makedirs(out_dir, exist_ok=True)
    out_name = f"{sample}_{assembler}_{config}_{lineage}"
    busco_subdir = os.path.join(out_dir, out_name)

    # If there's an existing short_summary.*.json, use it
    json_candidates = glob.glob(os.path.join(busco_subdir, "short_summary.*.json"))
    if json_candidates and not force:
        logger.info(f"✅ Resuming: found existing BUSCO {lineage} in {busco_subdir}")
        json_file = json_candidates[0]
    else:
        cmd = (
            f"busco -i {fasta_file} -o {out_name} -l {lineage} -m genome "
            f"--cpu {threads} --out_path {out_dir} --force"
        )
        ok, _ = run_in_env(env_name, cmd, desc=f"BUSCO {lineage} for {sample}/{assembler}/{config}")
        if not ok:
            if state_dir:
                update_state(state_dir, sample, assembler, config, step, "failed", error=f"busco_{lineage}_failed")
            return None

        json_candidates = glob.glob(os.path.join(busco_subdir, "short_summary.*.json"))
        if not json_candidates:
            logger.error(f"No BUSCO JSON found for lineage={lineage} in {busco_subdir}")
            if state_dir:
                update_state(state_dir, sample, assembler, config, step, "failed", error=f"busco_{lineage}_missing_json")
            return None
        json_file = json_candidates[0]

    busco_results = parse_busco_json(json_file)
    if state_dir:
        update_state(state_dir, sample, assembler, config, step, "completed",
                     path=json_file,
                     result=busco_results,
                     timestamp=datetime.now().isoformat())
    return json_file

# --------------------------------------------------------------------------------
# 8. Summaries
# --------------------------------------------------------------------------------
def create_comparison_summary(results, out_dir):
    """
    Create a DataFrame from 'results', sort by busco_mollusca DESC, busco_eukaryota DESC, n50 DESC.
    Save to CSV and MD.
    """
    df = pd.DataFrame(results)
    if df.empty:
        logger.warning("No results to summarize.")
        return

    df = df.sort_values(by=["busco_mollusca", "busco_eukaryota", "n50"], ascending=False)

    best = df.iloc[0]
    csv_path = os.path.join(out_dir, "assembly_comparison_summary.csv")
    df.to_csv(csv_path, index=False)

    md_path = os.path.join(out_dir, "assembly_comparison_summary.md")
    with open(md_path, 'w') as f:
        f.write("# Assembly Comparison Summary\n\n")
        f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write("Sorted by **Mollusca completeness** desc, **Eukaryota completeness** desc, then **N50** desc.\n\n")
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
    Collect results from the pipeline state to create final summary.
    Returns a list of dicts with keys:
      sample, assembler, config, busco_mollusca, missing_mollusca, busco_eukaryota,
      missing_eukaryota, n50, num_contigs, total_length, gc
    """
    st = read_state(output_dir)
    results = []
    for sample, sample_data in st.items():
        if "assemblers" not in sample_data:
            continue
        for asm_key, asm_data in sample_data["assemblers"].items():
            parts = asm_key.split("_", 1)
            if len(parts) < 2:
                continue
            assembler, config = parts

            assembly_data = asm_data.get("assembly", {})
            if assembly_data.get("status") != "completed":
                continue

            busco_moll = asm_data.get("busco_mollusca", {})
            busco_euk  = asm_data.get("busco_eukaryota", {})
            gc_data    = asm_data.get("gc", {})

            row = {
                "sample": sample,
                "assembler": assembler,
                "config": config,
                "busco_mollusca": 0.0,
                "missing_mollusca": 0.0,
                "busco_eukaryota": 0.0,
                "missing_eukaryota": 0.0,
                "n50": 0,
                "num_contigs": 0,
                "total_length": 0,
                "gc": 0.0
            }
            if busco_moll.get("status") == "completed":
                r = busco_moll.get("result", {})
                row["busco_mollusca"] = r.get("complete_percent", 0.0)
                row["missing_mollusca"] = r.get("missing_percent", 0.0)
                if row["n50"] == 0:
                    row["n50"] = r.get("n50", 0)
                if row["num_contigs"] == 0:
                    row["num_contigs"] = r.get("num_contigs", 0)
                if row["total_length"] == 0:
                    row["total_length"] = r.get("total_length", 0)
            if busco_euk.get("status") == "completed":
                r = busco_euk.get("result", {})
                row["busco_eukaryota"] = r.get("complete_percent", 0.0)
                row["missing_eukaryota"] = r.get("missing_percent", 0.0)
                if row["n50"] == 0:
                    row["n50"] = r.get("n50", 0)
                if row["num_contigs"] == 0:
                    row["num_contigs"] = r.get("num_contigs", 0)
                if row["total_length"] == 0:
                    row["total_length"] = r.get("total_length", 0)
            if gc_data.get("status") == "completed":
                row["gc"] = gc_data.get("value", 0.0)

            results.append(row)
    return results

# --------------------------------------------------------------------------------
# 9. Main
# --------------------------------------------------------------------------------
def main():
    # 1) Create conda envs with mamba if missing
    setup_conda_envs()

    # 2) Argparse
    parser = argparse.ArgumentParser(
        description="All-vs-All pipeline for SOAP(Multi-step+GapCloser), MaSuRCA, Captus, with BUSCO and GC%. With resume functionality."
    )
    parser.add_argument("-i", "--input_dir", required=True, help="Directory with cleaned FASTQ (R1/R2) files")
    parser.add_argument("-o", "--output_dir", default="Assembly_Comparison", help="Output directory")
    parser.add_argument("-t", "--threads", type=int, default=16, help="Number of CPU threads")
    parser.add_argument("--clean-after", action="store_true", help="Delete intermediate build folders after completion")
    parser.add_argument("--force", action="store_true", help="Force re-run all steps, ignoring previous results")
    parser.add_argument("--force-assembly", action="store_true", help="Force re-run assembly steps only")
    parser.add_argument("--force-busco", action="store_true", help="Force re-run BUSCO steps only")
    parser.add_argument("--force-gc", action="store_true", help="Force re-run GC calculation only")
    parser.add_argument("--show-state", action="store_true", help="Show current state and exit")
    parser.add_argument("--reset-state", action="store_true", help="Reset pipeline state before running")
    args = parser.parse_args()

    log_fp = setup_logging(args.output_dir)
    logger.info("=== Assembly Pipeline: SOAP + MaSuRCA + Captus ===")
    logger.info(f"Input dir: {args.input_dir}")
    logger.info(f"Output dir: {args.output_dir}")
    logger.info(f"Threads: {args.threads}")
    logger.info(f"Logs -> {log_fp}")

    # Reset state if requested
    if args.reset_state:
        sf = state_file_path(args.output_dir)
        if os.path.exists(sf):
            os.unlink(sf)
            logger.info(f"🧹 Reset pipeline state file: {sf}")

    # Show state if requested
    if args.show_state:
        st = read_state(args.output_dir)
        if not st:
            logger.info("No pipeline state found or state is empty.")
        else:
            logger.info("Current pipeline state:")
            for sample, data in st.items():
                logger.info(f"Sample: {sample}")
                for asm_key, asm_data in data.get("assemblers", {}).items():
                    logger.info(f"  {asm_key}:")
                    for step, step_data in asm_data.items():
                        status = step_data.get("status", "unknown")
                        status_symbol = "✅" if status == "completed" else "❌"
                        logger.info(f"    {step}: {status_symbol} {status}")
            # Summaries
            existing_results = collect_results_from_state(args.output_dir)
            if existing_results:
                df = pd.DataFrame(existing_results)
                df = df.sort_values(by=["busco_mollusca", "busco_eukaryota", "n50"], ascending=False)
                logger.info("\nResults summary (top 5):")
                pd.set_option('display.max_columns', None)
                logger.info(df.head(5).to_string())
        sys.exit(0)

    # Make subdirs
    asm_dir   = os.path.join(args.output_dir, "assemblies")
    busco_dir = os.path.join(args.output_dir, "busco")
    os.makedirs(asm_dir, exist_ok=True)
    os.makedirs(busco_dir, exist_ok=True)

    # Gather pairs
    r1_files = sorted(glob.glob(os.path.join(args.input_dir, "*_R1*.f*q*")))
    samples = {}
    for r1 in r1_files:
        base = os.path.basename(r1)
        sample_name = base.split("_R1")[0]
        r2_candidates = glob.glob(os.path.join(args.input_dir, sample_name + "_R2*.f*q*"))
        if r2_candidates:
            samples[sample_name] = (r1, r2_candidates[0])
        else:
            logger.warning(f"No R2 found for {r1}")

    if not samples:
        logger.error("No paired samples found. Exiting.")
        sys.exit(1)

    # We store final results in a list of dicts
    all_results = []
    lineages = ["mollusca_odb10", "eukaryota_odb10"]

    # Check if we are resuming
    state_exists = os.path.exists(state_file_path(args.output_dir)) and not args.force and not args.reset_state
    if state_exists:
        logger.info("🔄 Resuming pipeline from previous run.")
        existing_results = collect_results_from_state(args.output_dir)
        all_results.extend(existing_results)
        completed_configs = set(
            f"{r['sample']}_{r['assembler']}_{r['config']}" for r in existing_results
        )
        logger.info(f"Found {len(existing_results)} completed results in state.")
    else:
        completed_configs = set()

    # Determine forced re-runs
    force_assembly = args.force or args.force_assembly
    force_busco    = args.force or args.force_busco
    force_gc       = args.force or args.force_gc

    # Pipeline
    for sample_name, (r1, r2) in samples.items():
        logger.info(f"\n=== Sample: {sample_name} ===")

        # 1) SOAP k=63
        ckey = f"{sample_name}_SOAP_k63"
        if ckey not in completed_configs or force_assembly:
            soap63 = run_soapdenovo2_with_gapcloser(r1, r2, asm_dir, sample_name, 63,
                                                    threads=args.threads, state_dir=args.output_dir,
                                                    force=force_assembly)
            if soap63:
                gc_val = get_gc_content("soap_env", soap63, args.output_dir, sample_name, "SOAP", "k63", force=force_gc)
                busco_data = {}
                for lin in lineages:
                    bdir = os.path.join(busco_dir, f"soap_k63_{lin}")
                    busco_json = run_busco("soap_env", soap63, bdir, lin, args.threads, sample_name, "SOAP", "k63",
                                           state_dir=args.output_dir, force=force_busco)
                    if busco_json:
                        busco_data[lin] = parse_busco_json(busco_json)
                # Collect
                moll = busco_data.get("mollusca_odb10", {})
                euk  = busco_data.get("eukaryota_odb10", {})
                row = {
                    "sample": sample_name,
                    "assembler": "SOAP",
                    "config": "k63",
                    "busco_mollusca": moll.get("complete_percent", 0.0),
                    "missing_mollusca": moll.get("missing_percent", 0.0),
                    "busco_eukaryota": euk.get("complete_percent", 0.0),
                    "missing_eukaryota": euk.get("missing_percent", 0.0),
                    "n50": max(moll.get("n50",0), euk.get("n50",0)),
                    "num_contigs": max(moll.get("num_contigs",0), euk.get("num_contigs",0)),
                    "total_length": max(moll.get("total_length",0), euk.get("total_length",0)),
                    "gc": gc_val
                }
                # remove any old row for this config
                all_results = [r for r in all_results if not(r["sample"]==sample_name and r["assembler"]=="SOAP" and r["config"]=="k63")]
                all_results.append(row)
        else:
            logger.info(f"✅ Skipping SOAP k=63 for {sample_name}; already completed.")

        # 2) SOAP k=127
        ckey = f"{sample_name}_SOAP_k127"
        if ckey not in completed_configs or force_assembly:
            soap127 = run_soapdenovo2_with_gapcloser(r1, r2, asm_dir, sample_name, 127,
                                                     threads=args.threads, state_dir=args.output_dir,
                                                     force=force_assembly)
            if soap127:
                gc_val = get_gc_content("soap_env", soap127, args.output_dir, sample_name, "SOAP", "k127", force=force_gc)
                busco_data = {}
                for lin in lineages:
                    bdir = os.path.join(busco_dir, f"soap_k127_{lin}")
                    busco_json = run_busco("soap_env", soap127, bdir, lin, args.threads, sample_name, "SOAP", "k127",
                                           state_dir=args.output_dir, force=force_busco)
                    if busco_json:
                        busco_data[lin] = parse_busco_json(busco_json)
                moll = busco_data.get("mollusca_odb10", {})
                euk  = busco_data.get("eukaryota_odb10", {})
                row = {
                    "sample": sample_name,
                    "assembler": "SOAP",
                    "config": "k127",
                    "busco_mollusca": moll.get("complete_percent", 0.0),
                    "missing_mollusca": moll.get("missing_percent", 0.0),
                    "busco_eukaryota": euk.get("complete_percent", 0.0),
                    "missing_eukaryota": euk.get("missing_percent", 0.0),
                    "n50": max(moll.get("n50",0), euk.get("n50",0)),
                    "num_contigs": max(moll.get("num_contigs",0), euk.get("num_contigs",0)),
                    "total_length": max(moll.get("total_length",0), euk.get("total_length",0)),
                    "gc": gc_val
                }
                all_results = [r for r in all_results if not(r["sample"]==sample_name and r["assembler"]=="SOAP" and r["config"]=="k127")]
                all_results.append(row)
        else:
            logger.info(f"✅ Skipping SOAP k=127 for {sample_name}; already completed.")

        # 3) MaSuRCA "permissive"
        ckey = f"{sample_name}_MaSuRCA_permissive"
        if ckey not in completed_configs or force_assembly:
            ms_perm = run_masurca(r1, r2, asm_dir, sample_name, "permissive",
                                  threads=args.threads, state_dir=args.output_dir, force=force_assembly)
            if ms_perm:
                gc_val = get_gc_content("masurca_env", ms_perm, args.output_dir, sample_name, "MaSuRCA", "permissive", force=force_gc)
                busco_data = {}
                for lin in lineages:
                    bdir = os.path.join(busco_dir, f"masurca_permissive_{lin}")
                    busco_json = run_busco("masurca_env", ms_perm, bdir, lin, args.threads, sample_name, "MaSuRCA", "permissive",
                                           state_dir=args.output_dir, force=force_busco)
                    if busco_json:
                        busco_data[lin] = parse_busco_json(busco_json)
                moll = busco_data.get("mollusca_odb10", {})
                euk  = busco_data.get("eukaryota_odb10", {})
                row = {
                    "sample": sample_name,
                    "assembler": "MaSuRCA",
                    "config": "permissive",
                    "busco_mollusca": moll.get("complete_percent", 0.0),
                    "missing_mollusca": moll.get("missing_percent", 0.0),
                    "busco_eukaryota": euk.get("complete_percent", 0.0),
                    "missing_eukaryota": euk.get("missing_percent", 0.0),
                    "n50": max(moll.get("n50",0), euk.get("n50",0)),
                    "num_contigs": max(moll.get("num_contigs",0), euk.get("num_contigs",0)),
                    "total_length": max(moll.get("total_length",0), euk.get("total_length",0)),
                    "gc": gc_val
                }
                all_results = [r for r in all_results if not(r["sample"]==sample_name and r["assembler"]=="MaSuRCA" and r["config"]=="permissive")]
                all_results.append(row)
        else:
            logger.info(f"✅ Skipping MaSuRCA permissive for {sample_name}; already completed.")

        # 4) MaSuRCA "stringent"
        ckey = f"{sample_name}_MaSuRCA_stringent"
        if ckey not in completed_configs or force_assembly:
            ms_str = run_masurca(r1, r2, asm_dir, sample_name, "stringent",
                                 threads=args.threads, state_dir=args.output_dir, force=force_assembly)
            if ms_str:
                gc_val = get_gc_content("masurca_env", ms_str, args.output_dir, sample_name, "MaSuRCA", "stringent", force=force_gc)
                busco_data = {}
                for lin in lineages:
                    bdir = os.path.join(busco_dir, f"masurca_stringent_{lin}")
                    busco_json = run_busco("masurca_env", ms_str, bdir, lin, args.threads, sample_name, "MaSuRCA", "stringent",
                                           state_dir=args.output_dir, force=force_busco)
                    if busco_json:
                        busco_data[lin] = parse_busco_json(busco_json)
                moll = busco_data.get("mollusca_odb10", {})
                euk  = busco_data.get("eukaryota_odb10", {})
                row = {
                    "sample": sample_name,
                    "assembler": "MaSuRCA",
                    "config": "stringent",
                    "busco_mollusca": moll.get("complete_percent", 0.0),
                    "missing_mollusca": moll.get("missing_percent", 0.0),
                    "busco_eukaryota": euk.get("complete_percent", 0.0),
                    "missing_eukaryota": euk.get("missing_percent", 0.0),
                    "n50": max(moll.get("n50",0), euk.get("n50",0)),
                    "num_contigs": max(moll.get("num_contigs",0), euk.get("num_contigs",0)),
                    "total_length": max(moll.get("total_length",0), euk.get("total_length",0)),
                    "gc": gc_val
                }
                all_results = [r for r in all_results if not(r["sample"]==sample_name and r["assembler"]=="MaSuRCA" and r["config"]=="stringent")]
                all_results.append(row)
        else:
            logger.info(f"✅ Skipping MaSuRCA stringent for {sample_name}; already completed.")

        # 5) Captus
        ckey = f"{sample_name}_Captus_WGS"
        if ckey not in completed_configs or force_assembly:
            captus_fa = run_captus_assembly(r1, r2, asm_dir, sample_name, 
                                            threads=args.threads, state_dir=args.output_dir, force=force_assembly)
            if captus_fa:
                gc_val = get_gc_content("captus_env", captus_fa, args.output_dir, sample_name, "Captus", "WGS", force=force_gc)
                busco_data = {}
                for lin in lineages:
                    bdir = os.path.join(busco_dir, f"captus_wgs_{lin}")
                    busco_json = run_busco("captus_env", captus_fa, bdir, lin, args.threads, sample_name, "Captus", "WGS",
                                           state_dir=args.output_dir, force=force_busco)
                    if busco_json:
                        busco_data[lin] = parse_busco_json(busco_json)
                moll = busco_data.get("mollusca_odb10", {})
                euk  = busco_data.get("eukaryota_odb10", {})
                row = {
                    "sample": sample_name,
                    "assembler": "Captus",
                    "config": "WGS",
                    "busco_mollusca": moll.get("complete_percent", 0.0),
                    "missing_mollusca": moll.get("missing_percent", 0.0),
                    "busco_eukaryota": euk.get("complete_percent", 0.0),
                    "missing_eukaryota": euk.get("missing_percent", 0.0),
                    "n50": max(moll.get("n50",0), euk.get("n50",0)),
                    "num_contigs": max(moll.get("num_contigs",0), euk.get("num_contigs",0)),
                    "total_length": max(moll.get("total_length",0), euk.get("total_length",0)),
                    "gc": gc_val
                }
                all_results = [r for r in all_results if not(r["sample"]==sample_name and r["assembler"]=="Captus" and r["config"]=="WGS")]
                all_results.append(row)
        else:
            logger.info(f"✅ Skipping Captus WGS for {sample_name}; already completed.")

    # Summaries
    create_comparison_summary(all_results, args.output_dir)

    if args.clean_after:
        for sample_name in samples:
            for d in [
                os.path.join(asm_dir, f"{sample_name}_soap_k63"),
                os.path.join(asm_dir, f"{sample_name}_soap_k127"),
                os.path.join(asm_dir, f"{sample_name}_masurca_permissive"),
                os.path.join(asm_dir, f"{sample_name}_masurca_stringent"),
                os.path.join(asm_dir, f"{sample_name}_captus_wgs_run")
            ]:
                if os.path.exists(d):
                    shutil.rmtree(d)
                    logger.info(f"🧹 Cleaned up temp folder: {d}")

    logger.info("Pipeline complete.")

if __name__ == "__main__":
    main()
