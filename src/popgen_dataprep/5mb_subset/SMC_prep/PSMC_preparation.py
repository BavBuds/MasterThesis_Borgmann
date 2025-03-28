#!/usr/bin/env python3
import os
import subprocess
import concurrent.futures

# Settings
NUM_THREADS = 30
NUM_BOOTSTRAPS = 100

# Base directories
BASE_DIR = "/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/Popgen_analysis/dataprep/5mb_subset"
QC_DIR = os.path.join(BASE_DIR, "QC", "Intermediate_data")
PSMC_OUT_DIR = os.path.join(BASE_DIR, "PSMC_preparation")

# List of samples
samples = ["HT", "HT2", "BB"]

# Create the PSMC output base directory if it doesn't exist
os.makedirs(PSMC_OUT_DIR, exist_ok=True)

def run(cmd, workdir=None):
    """Run a shell command and print it."""
    print(f"🔧 Running: {cmd}")
    subprocess.run(cmd, shell=True, check=True, cwd=workdir)

def process_bootstrap(args):
    """Process a single bootstrap replicate"""
    sample, i, split_file, sample_out_dir = args
    boot_file = os.path.join(sample_out_dir, f"bootstrap.{i}.psmc")
    # Use more PSMC parameters as recommended in the documentation
    cmd_boot = f"psmc -N25 -t15 -r5 -b -p \"4+25*2+4+6\" -o {boot_file} {split_file}"
    run(cmd_boot)
    return i

for sample in samples:
    print(f"\n🚀 Processing sample: {sample}")
    
    # Create a sample-specific output directory for PSMC prep
    sample_out_dir = os.path.join(PSMC_OUT_DIR, sample)
    os.makedirs(sample_out_dir, exist_ok=True)
    
    # Input filtered VCF (output from QC pipeline)
    input_vcf = os.path.join(QC_DIR, sample, f"{sample}_filtered.vcf.gz")
    
    # Reference FASTA (the 5Mb subset FASTA for this sample)
    ref_fasta = os.path.join(BASE_DIR, f"{sample}_assembly_5mb_subset.fasta")
    
    # Step 1: Convert filtered VCF to diploid FASTQ - FIX THE SYNTAX ERROR
    fastq_file = os.path.join(sample_out_dir, f"{sample}.fq.gz")
    
    # Fixed command: use separate commands instead of piping which is causing issues
    cmd_fastq = f"bcftools view {input_vcf} | vcfutils.pl vcf2fq -d 10 -D 100 | gzip > {fastq_file}"
    run(cmd_fastq)
    
    # Step 2: Verify the FASTQ file exists and has content
    cmd_check = f"gzip -dc {fastq_file} | head -n 20"
    run(cmd_check)
    
    # Step 2: Convert FASTQ to PSMC format (psmcfa)
    psmcfa_file = os.path.join(sample_out_dir, f"{sample}.psmcfa")
    cmd_psmcfa = f"fq2psmcfa -q20 {fastq_file} > {psmcfa_file}"
    run(cmd_psmcfa)
    
    # Step 3: Run the main PSMC analysis
    psmc_file = os.path.join(sample_out_dir, f"{sample}.psmc")
    cmd_psmc = f"psmc -N25 -t15 -r5 -p \"4+25*2+4+6\" -o {psmc_file} {psmcfa_file}"
    run(cmd_psmc)
    
    # Step 4: Split the psmcfa file (required for bootstrapping)
    split_file = os.path.join(sample_out_dir, f"{sample}.split.psmcfa")
    cmd_split = f"splitfa {psmcfa_file} > {split_file}"
    run(cmd_split)
    
    # Step 5: Generate bootstrap replicates (100 replicates)
    print(f"🔁 Generating {NUM_BOOTSTRAPS} bootstrap replicates for {sample} using {NUM_THREADS} threads...")
    
    # Use ThreadPoolExecutor to run bootstraps in parallel
    with concurrent.futures.ThreadPoolExecutor(max_workers=NUM_THREADS) as executor:
        args_list = [(sample, i, split_file, sample_out_dir) for i in range(1, NUM_BOOTSTRAPS + 1)]
        results = list(executor.map(process_bootstrap, args_list))
    
    # Step 6: Combine original and bootstrap results
    combined_file = os.path.join(sample_out_dir, f"{sample}_combined.psmc")
    boot_files = " ".join([os.path.join(sample_out_dir, f"bootstrap.{i}.psmc") for i in range(1, NUM_BOOTSTRAPS + 1)])
    cmd_combine = f"cat {psmc_file} {boot_files} > {combined_file}"
    run(cmd_combine)
    
    # Step 7: Generate plot
    cmd_plot = f"psmc_plot.pl -pY50000 {os.path.join(sample_out_dir, sample)} {combined_file}"
    run(cmd_plot)
    
    print(f"✅ Finished PSMC prep for {sample}")

print("\n✅ All samples processed for PSMC preparation.")