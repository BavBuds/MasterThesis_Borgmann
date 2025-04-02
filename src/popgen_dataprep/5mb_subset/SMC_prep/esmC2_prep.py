#!/usr/bin/env python3
"""
eSMC2 Preparation Script

This script processes filtered BAM files (output from the QC pipeline) to create input files
for eSMC2 analysis. It automatically adds read group information if missing from BAM files.

The workflow includes:
1. Check and fix BAM read groups
2. Automatic reference indexing (.fai and .dict creation)
3. HaplotypeCaller to generate GVCFs
4. CombineGVCFs for specified cohorts
5. GenotypeGVCFs on combined GVCFs
6. Creating SNP VCFs and callable sites VCFs
7. Generating mask files
8. Creating Multihetsep files for eSMC2

Author: Max Borgmann
Date: April 2025
"""

import os
import subprocess
import sys
import argparse
import shutil
import glob
import re
import time
import tempfile
from datetime import datetime

def check_dependencies():
    """Check if required software is installed."""
    required_tools = {
        "gatk": "GATK",
        "bcftools": "Bcftools",
        "vcf2bed": "Bedops",
        "bedtools": "Bedtools",
        "samtools": "Samtools",
        "wget": "Wget"
    }
    
    missing_tools = []
    for cmd, name in required_tools.items():
        if shutil.which(cmd) is None:
            missing_tools.append(name)
    
    if missing_tools:
        print(f"⚠️  Missing required tools: {', '.join(missing_tools)}")
        print("Please install them before running this script.")
        return False
    
    print("✅ All required tools are installed.")
    return True

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Prepare input files for eSMC2 analysis from filtered BAM files.")
    parser.add_argument("--input-fastas-dir", required=True, 
                        help="Directory containing sample-specific reference fasta files (Input_fastas directory)")
    parser.add_argument("--qc-dir", required=True, 
                        help="Directory containing filtered BAM files from QC pipeline")
    parser.add_argument("--output-dir", required=True, 
                        help="Output directory for eSMC2 preparation")
    parser.add_argument("--depth-filter", type=int, default=30, 
                        help="Minimum depth filter (default: 30)")
    parser.add_argument("--memory", default="4g", 
                        help="Memory allocation for Java (default: 4g)")
    parser.add_argument("--threads", type=int, default=8, 
                        help="Number of threads (default: 8)")
    parser.add_argument("--region", 
                        help="Specific genomic region to analyze (optional)")
    parser.add_argument("--create-subset", action="store_true", 
                        help="Create a subset of the Multihetsep file")
    parser.add_argument("--subset-start", type=int, default=40000000, 
                        help="Start position for subset (default: 40000000)")
    parser.add_argument("--subset-end", type=int, default=45000000, 
                        help="End position for subset (default: 45000000)")
    parser.add_argument("--log-file", default="esmc2_prep.log", 
                        help="Log file name (default: esmc2_prep.log)")
    return parser.parse_args()

def setup_logging(log_file):
    """Set up logging to both console and file."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] eSMC2 preparation pipeline started")
    print(f"[{timestamp}] Log will be written to: {log_file}")
    
    # Redirect stdout and stderr to both console and file
    class Logger:
        def __init__(self, filename):
            self.terminal = sys.stdout
            self.log = open(filename, "w")
            
        def write(self, message):
            self.terminal.write(message)
            self.log.write(message)
            self.flush()
            
        def flush(self):
            self.terminal.flush()
            self.log.flush()
    
    sys.stdout = Logger(log_file)
    sys.stderr = sys.stdout
    
    return timestamp

def log(message, level="INFO"):
    """Log a message with timestamp."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {level}: {message}")

def run(cmd, workdir=None):
    """Execute a shell command and log the output."""
    log(f"Running: {cmd}", "CMD")
    try:
        start_time = time.time()
        process = subprocess.run(
            cmd, 
            shell=True, 
            check=True, 
            cwd=workdir,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True
        )
        elapsed_time = time.time() - start_time
        log(f"Command completed in {elapsed_time:.2f} seconds", "INFO")
        if process.stdout.strip():
            log(f"stdout: {process.stdout.strip()}", "DEBUG")
        if process.stderr.strip():
            log(f"stderr: {process.stderr.strip()}", "DEBUG")
        return True, process.stdout
    except subprocess.CalledProcessError as e:
        log(f"Command failed with exit code {e.returncode}", "ERROR")
        log(f"stdout: {e.stdout.strip() if e.stdout else 'None'}", "ERROR")
        log(f"stderr: {e.stderr.strip() if e.stderr else 'None'}", "ERROR")
        return False, None

def index_reference_fasta(fasta_path):
    """
    Create necessary index files for a reference FASTA:
    1. .fai index using samtools faidx
    2. .dict dictionary using GATK CreateSequenceDictionary
    
    Returns True if indexing was successful or indices already exist
    """
    fai_path = fasta_path + ".fai"
    dict_path = os.path.splitext(fasta_path)[0] + ".dict"
    
    # Check if indexes already exist
    fai_exists = os.path.exists(fai_path)
    dict_exists = os.path.exists(dict_path)
    
    if fai_exists and dict_exists:
        log(f"Reference indices already exist for {fasta_path}", "INFO")
        return True
    
    # Create .fai index if needed
    if not fai_exists:
        log(f"Creating .fai index for {fasta_path}", "INFO")
        faidx_cmd = f"samtools faidx {fasta_path}"
        if not run(faidx_cmd)[0]:
            log(f"Failed to create .fai index for {fasta_path}", "ERROR")
            return False
    
    # Create .dict dictionary if needed
    if not dict_exists:
        log(f"Creating .dict dictionary for {fasta_path}", "INFO")
        dict_cmd = f"gatk CreateSequenceDictionary -R {fasta_path}"
        if not run(dict_cmd)[0]:
            log(f"Failed to create .dict dictionary for {fasta_path}", "ERROR")
            return False
    
    return True

def check_and_fix_read_groups(bam_file, sample_name):
    """
    Check if BAM file has read group tags with sample name.
    If not, add them and create a new BAM file.
    
    Returns the path to the BAM file with proper read groups.
    """
    # Check if BAM has read groups with SM tag
    cmd = f"samtools view -H {bam_file} | grep '^@RG' | grep -o 'SM:[^\\t]*' | head -1"
    success, output = run(cmd)
    
    if success and output and 'SM:' in output:
        # Extract existing sample name from the SM tag
        existing_sample = output.strip().split('SM:')[1]
        log(f"BAM file already has read group with sample name: {existing_sample}", "INFO")
        return bam_file, existing_sample
    
    # BAM doesn't have read groups or SM tag, add them
    log(f"BAM file doesn't have proper read groups. Adding @RG with sample name: {sample_name}", "INFO")
    
    # Create output directory
    output_dir = os.path.join(os.path.dirname(bam_file), "fixed_bams")
    os.makedirs(output_dir, exist_ok=True)
    
    # Define output BAM with read groups
    fixed_bam = os.path.join(output_dir, os.path.basename(bam_file).replace(".bam", ".with_rg.bam"))
    
    # Skip if fixed BAM already exists
    if os.path.exists(fixed_bam):
        log(f"Fixed BAM already exists: {fixed_bam}", "INFO")
        return fixed_bam, sample_name
    
    # Add read groups
    rg_cmd = (
        f"samtools addreplacerg -r '@RG\\tID:{sample_name}\\tSM:{sample_name}\\tPL:ILLUMINA' "
        f"-o {fixed_bam} {bam_file}"
    )
    
    if run(rg_cmd)[0]:
        log(f"Successfully added read groups to BAM: {fixed_bam}", "INFO")
        
        # Index the fixed BAM
        index_cmd = f"samtools index {fixed_bam}"
        if run(index_cmd)[0]:
            log(f"Indexed fixed BAM: {fixed_bam}", "INFO")
            return fixed_bam, sample_name
        else:
            log(f"Failed to index fixed BAM: {fixed_bam}", "ERROR")
            return bam_file, sample_name
    else:
        log(f"Failed to add read groups to BAM: {bam_file}", "ERROR")
        return bam_file, sample_name

def find_reference_fastas(input_fastas_dir):
    """Find sample-specific reference fasta files and index them."""
    reference_fastas = {}
    
    # Check for each sample's subdirectory
    for sample_dir in glob.glob(os.path.join(input_fastas_dir, "*")):
        if os.path.isdir(sample_dir):
            sample_name = os.path.basename(sample_dir)
            
            # Look for sample_assembly_5mb_subset.fasta
            fasta_path = os.path.join(sample_dir, f"{sample_name}_assembly_5mb_subset.fasta")
            if os.path.exists(fasta_path):
                log(f"Found reference fasta for {sample_name}: {fasta_path}")
                
                # Index the reference fasta
                if index_reference_fasta(fasta_path):
                    reference_fastas[sample_name] = fasta_path
                else:
                    log(f"Skipping {sample_name} due to indexing failure", "WARNING")
            else:
                # Also try looking for other fasta files
                fasta_files = glob.glob(os.path.join(sample_dir, "*.fasta")) + glob.glob(os.path.join(sample_dir, "*.fa"))
                if fasta_files:
                    fasta_path = fasta_files[0]
                    log(f"Found alternative reference fasta for {sample_name}: {fasta_path}")
                    
                    # Index the reference fasta
                    if index_reference_fasta(fasta_path):
                        reference_fastas[sample_name] = fasta_path
                    else:
                        log(f"Skipping {sample_name} due to indexing failure", "WARNING")
                else:
                    log(f"No reference fasta found for {sample_name}", "WARNING")
    
    return reference_fastas

def discover_filtered_bams(qc_dir):
    """Discover filtered BAM files in the QC output directory."""
    filtered_bams = {}
    
    # Looking in QC/Intermediate_data/{sample}/{sample}_filtered.bam
    intermediate_dir = os.path.join(qc_dir, "Intermediate_data")
    if not os.path.exists(intermediate_dir):
        log(f"Intermediate data directory not found: {intermediate_dir}", "ERROR")
        return filtered_bams
    
    # Check each sample subdirectory
    for sample_dir in glob.glob(os.path.join(intermediate_dir, "*")):
        if os.path.isdir(sample_dir):
            sample_name = os.path.basename(sample_dir)
            # Look for sample_filtered.bam
            filtered_bam = os.path.join(sample_dir, f"{sample_name}_filtered.bam")
            # Also check for sample_sorted.bam as alternative
            sorted_bam = os.path.join(sample_dir, f"{sample_name}_sorted.bam")
            
            if os.path.exists(filtered_bam):
                filtered_bams[sample_name] = filtered_bam
                log(f"Found filtered BAM for {sample_name}: {filtered_bam}")
            elif os.path.exists(sorted_bam):
                filtered_bams[sample_name] = sorted_bam
                log(f"Found sorted BAM for {sample_name}: {sorted_bam}")
    
    return filtered_bams

def prompt_user_for_cohorts(available_samples, reference_fastas):
    """Prompt user to define cohorts from available samples."""
    print("\n📋 Available samples:")
    for i, sample in enumerate(sorted(available_samples), 1):
        ref_status = "✅" if sample in reference_fastas else "❌"
        print(f"  {i}. {sample} {ref_status}")
    
    if len(reference_fastas) < len(available_samples):
        print("\n⚠️  Warning: Some samples don't have reference fastas and will be skipped.")
    
    cohorts = {}
    while True:
        cohort_name = input("\n🔍 Enter a name for a new cohort (or press Enter to finish): ").strip()
        if not cohort_name:
            break
        
        print(f"\n🔍 Select samples for cohort '{cohort_name}':")
        print("Enter sample numbers separated by commas (e.g., 1,3), or 'all' for all samples")
        selection = input("Selection: ").strip()
        
        selected_samples = []
        if selection.lower() == 'all':
            selected_samples = list(available_samples)
        else:
            try:
                sorted_samples = sorted(available_samples)
                indices = [int(idx.strip()) - 1 for idx in selection.split(',')]
                selected_samples = [sorted_samples[i] for i in indices if 0 <= i < len(sorted_samples)]
            except ValueError:
                print("⚠️ Invalid selection. Please enter numbers separated by commas.")
                continue
        
        # Filter out samples without reference fastas
        valid_samples = [s for s in selected_samples if s in reference_fastas]
        if not valid_samples:
            print("⚠️ None of the selected samples have reference fastas. Please try again.")
            continue
        
        if len(valid_samples) < len(selected_samples):
            print(f"⚠️ Warning: {len(selected_samples) - len(valid_samples)} samples were skipped due to missing reference fastas.")
        
        if valid_samples:
            # For a cohort, use the reference fasta of the first sample in the cohort
            reference_for_cohort = reference_fastas[valid_samples[0]]
            print(f"Using reference fasta from {valid_samples[0]} for cohort {cohort_name}: {reference_for_cohort}")
            
            cohorts[cohort_name] = {
                "samples": [(s, available_samples[s]) for s in valid_samples],
                "reference": reference_for_cohort
            }
            print(f"✅ Created cohort '{cohort_name}' with {len(valid_samples)} samples: {', '.join(valid_samples)}")
        else:
            print("⚠️ No valid samples selected for this cohort.")
    
    return cohorts

def main():
    """Main function to prepare input files for eSMC2 analysis."""
    if not check_dependencies():
        return 1
    
    args = parse_args()
    
    # Set up output directory and logging
    os.makedirs(args.output_dir, exist_ok=True)
    log_file = os.path.join(args.output_dir, args.log_file)
    start_time = setup_logging(log_file)
    
    # Print configuration
    log(f"Configuration:")
    log(f"  Input fastas directory: {args.input_fastas_dir}")
    log(f"  QC directory: {args.qc_dir}")
    log(f"  Output directory: {args.output_dir}")
    log(f"  Threads: {args.threads}")
    log(f"  Memory: {args.memory}")
    log(f"  Depth filter: {args.depth_filter}")
    if args.region:
        log(f"  Region: {args.region}")
    if args.create_subset:
        log(f"  Creating subset from {args.subset_start} to {args.subset_end}")
    
    # Ensure directories exist
    for dir_path in [args.input_fastas_dir, args.qc_dir]:
        if not os.path.isdir(dir_path):
            log(f"Directory does not exist: {dir_path}", "ERROR")
            return 1
    
    # Find and index sample-specific reference fastas
    log("Finding and indexing reference FASTA files...")
    reference_fastas = find_reference_fastas(args.input_fastas_dir)
    if not reference_fastas:
        log(f"No reference fasta files found in {args.input_fastas_dir}", "ERROR")
        return 1
    
    log(f"Found and indexed {len(reference_fastas)} reference fasta files")
    
    # Create GVCF output directory
    gvcf_dir = os.path.join(args.output_dir, "GVCFs")
    os.makedirs(gvcf_dir, exist_ok=True)
    
    # Find filtered BAM files from QC pipeline
    filtered_bams = discover_filtered_bams(args.qc_dir)
    if not filtered_bams:
        log("No filtered BAM files found in the QC directory", "ERROR")
        return 1
    
    log(f"Found {len(filtered_bams)} filtered BAM files")
    
    # Check and fix read groups in BAM files
    fixed_bams = {}
    sample_names = {}
    for sample_name, bam_file in filtered_bams.items():
        fixed_bam, actual_sample = check_and_fix_read_groups(bam_file, sample_name)
        fixed_bams[sample_name] = fixed_bam
        sample_names[sample_name] = actual_sample
    
    # Prompt user to define cohorts
    cohorts = prompt_user_for_cohorts(filtered_bams, reference_fastas)
    if not cohorts:
        log("No cohorts defined. Exiting.", "WARNING")
        return 1
    
    # Region parameter for GATK
    region_param = f"-L {args.region}" if args.region else ""
    region_name = args.region.replace(":", "_") if args.region else "wholeGenome"
    
    # Step 1: Run HaplotypeCaller on each BAM file
    gvcfs = {}
    for sample_name, bam_file in fixed_bams.items():
        if sample_name not in reference_fastas:
            log(f"Skipping HaplotypeCaller for {sample_name} - no reference fasta found", "WARNING")
            continue
            
        reference_fasta = reference_fastas[sample_name]
        actual_sample = sample_names[sample_name]
        log(f"Running HaplotypeCaller on sample: {sample_name} (using sample name: {actual_sample})")
        gvcf_out = os.path.join(gvcf_dir, f"{sample_name}.{region_name}.g.vcf.gz")
        
        # Skip if GVCF already exists
        if os.path.exists(gvcf_out):
            log(f"GVCF already exists for {sample_name}: {gvcf_out}")
            gvcfs[sample_name] = gvcf_out
            continue
        
        # Run HaplotypeCaller
        hc_cmd = (
            f'gatk --java-options "-Xmx{args.memory}" HaplotypeCaller '
            f'-R {reference_fasta} '
            f'-I {bam_file} '
            f'-O {gvcf_out} '
            f'-ERC BP_RESOLUTION '
            f'--output-mode EMIT_ALL_CONFIDENT_SITES '
            f'--native-pair-hmm-threads {args.threads} '
            f'{region_param}'
        )
        
        if run(hc_cmd)[0]:
            gvcfs[sample_name] = gvcf_out
            log(f"Successfully created GVCF for {sample_name}")
        else:
            log(f"Failed to create GVCF for {sample_name}", "ERROR")
    
    if not gvcfs:
        log("No GVCFs were created. Exiting.", "ERROR")
        return 1
    
    # Process each cohort
    for cohort_name, cohort_info in cohorts.items():
        log(f"Processing cohort: {cohort_name}")
        cohort_dir = os.path.join(args.output_dir, cohort_name)
        os.makedirs(cohort_dir, exist_ok=True)
        
        # Get the sample names, corresponding GVCF files, and reference fasta for this cohort
        sample_list = cohort_info["samples"]
        reference_fasta = cohort_info["reference"]
        cohort_samples = [(s, gvcfs.get(s)) for s, _ in sample_list if s in gvcfs]
        
        if not cohort_samples or None in [g for _, g in cohort_samples]:
            missing_samples = [s for s, _ in sample_list if s not in gvcfs]
            log(f"Missing GVCFs for samples in cohort {cohort_name}: {', '.join(missing_samples)}", "ERROR")
            log(f"Skipping cohort {cohort_name}", "ERROR")
            continue
        
        # 1. Combine GVCFs
        combined_gvcf = os.path.join(cohort_dir, f"{cohort_name}.{region_name}.g.vcf.gz")
        combine_cmd = f"gatk CombineGVCFs -R {reference_fasta}"
        for sample, gvcf in cohort_samples:
            combine_cmd += f" --variant {gvcf}"
        combine_cmd += f" -O {combined_gvcf}"
        if not run(combine_cmd)[0]:
            log(f"Failed to combine GVCFs for cohort {cohort_name}", "ERROR")
            continue
        
        # 2. Genotype GVCFs
        allsites_vcf = os.path.join(cohort_dir, f"{cohort_name}.allsites.geno.vcf.gz")
        genotype_cmd = (
            f'gatk --java-options "-Xmx{args.memory}" GenotypeGVCFs '
            f'-R {reference_fasta} '
            f'-V {combined_gvcf} '
            f'--include-non-variant-sites '
            f'-O {allsites_vcf}'
        )
        if not run(genotype_cmd)[0]:
            log(f"Failed to genotype GVCFs for cohort {cohort_name}", "ERROR")
            continue
        
        # 3. Create SNP-only VCF
        final_filtered_vcf = os.path.join(cohort_dir, f"{cohort_name}.final.filtered.vcf.gz")
        snp_view_cmd = (
            f'bcftools view {allsites_vcf} '
            f'--genotype ^miss '
            f'--apply-filters .,PASS '
            f'--include \'TYPE="snp"\' '
            f'-Oz -o {final_filtered_vcf}'
        )
        if not run(snp_view_cmd)[0]:
            log(f"Failed to create SNP-only VCF for cohort {cohort_name}", "ERROR")
            continue
        
        # Index the filtered VCF
        if not run(f"bcftools index --threads {args.threads} {final_filtered_vcf}")[0]:
            log(f"Failed to index filtered VCF for cohort {cohort_name}", "ERROR")
            continue
        
        # 4. Create callable sites VCF
        allsites_filtered_vcf = os.path.join(cohort_dir, f"{cohort_name}.final.allsites.geno.filtered.vcf.gz")
        mask_vcf_cmd = (
            f'bcftools view {allsites_vcf} '
            f'--genotype ^miss '
            f'--apply-filters .,PASS '
            f'--include \'TYPE="snp" && INFO/DP > {args.depth_filter} || TYPE="ref" && INFO/DP > {args.depth_filter}\' '
            f'-Oz -o {allsites_filtered_vcf}'
        )
        if not run(mask_vcf_cmd)[0]:
            log(f"Failed to create callable sites VCF for cohort {cohort_name}", "ERROR")
            continue
        
        # 5. Create mask file
        mask_bed = os.path.join(cohort_dir, f"{cohort_name}.final.mask.bed")
        vcf2bed_cmd = f"zcat {allsites_filtered_vcf} | vcf2bed > {mask_bed}"
        if not run(vcf2bed_cmd)[0]:
            log(f"Failed to convert VCF to BED for cohort {cohort_name}", "ERROR")
            continue
        
        merged_mask_bed = os.path.join(cohort_dir, f"{cohort_name}.final.mask.merged.bed")
        merge_cmd = f"bedtools merge -i {mask_bed} > {merged_mask_bed}"
        if not run(merge_cmd)[0]:
            log(f"Failed to merge BED intervals for cohort {cohort_name}", "ERROR")
            continue
        
        final_mask = merged_mask_bed + ".gz"
        if not run(f"gzip -f {merged_mask_bed}")[0]:
            log(f"Failed to compress mask BED for cohort {cohort_name}", "ERROR")
            continue
        
        # 6. Split VCF into individual samples
        sample_list_cmd = f"bcftools query -l {final_filtered_vcf}"
        success, output = run(sample_list_cmd)
        if success and output:
            sample_names_in_vcf = output.strip().split()
            log(f"Sample names in VCF: {', '.join(sample_names_in_vcf)}")
        else:
            log(f"Failed to extract sample names from {final_filtered_vcf}", "ERROR")
            continue
        
        per_sample_vcfs = []
        for s in sample_names_in_vcf:
            out_vcf = os.path.join(cohort_dir, f"{cohort_name}.{s}.vcf.gz")
            extract_cmd = f"bcftools view --samples {s} {final_filtered_vcf} --output-type z --output-file {out_vcf}"
            if run(extract_cmd)[0]:
                per_sample_vcfs.append(out_vcf)
                log(f"Created individual sample VCF for {s}")
            else:
                log(f"Failed to create individual sample VCF for {s}", "ERROR")
        
        if not per_sample_vcfs:
            log(f"No sample VCFs created for cohort {cohort_name}", "ERROR")
            continue
        
        # 7. Generate Multihetsep file
        gen_mhs_script = os.path.join(cohort_dir, "generate_multihetsep.py")
        if not os.path.exists(gen_mhs_script):
            log(f"Downloading generate_multihetsep.py script")
            wget_cmd = f"wget -O {gen_mhs_script} https://raw.githubusercontent.com/stschiff/msmc-tools/master/generate_multihetsep.py"
            if not run(wget_cmd)[0] or not run(f"chmod u+x {gen_mhs_script}")[0]:
                log(f"Failed to download or set permissions for generate_multihetsep.py", "ERROR")
                continue
        
        multihetsep_out = os.path.join(cohort_dir, f"{cohort_name}.{region_name}.mhs")
        mhs_cmd = f"python3 {gen_mhs_script} --mask={final_mask} " + " ".join(per_sample_vcfs) + f" > {multihetsep_out}"
        if not run(mhs_cmd)[0]:
            log(f"Failed to generate Multihetsep file for cohort {cohort_name}", "ERROR")
            continue
        
        # 8. Create subset if requested
        if args.create_subset:
            subset_out = os.path.join(
                cohort_dir, 
                f"{cohort_name}.{region_name}.{args.subset_start/1000000:.0f}Mbto{args.subset_end/1000000:.0f}Mb.subset.mhs"
            )
            subset_cmd = (
                f"awk -F'\\t' '$2>{args.subset_start} && $2<{args.subset_end}' {multihetsep_out} | "
                f"awk -v OFS='\\t' '{{$2 = $2 - {args.subset_start}; print}}' > {subset_out}"
            )
            if run(subset_cmd)[0]:
                log(f"Created subset file: {subset_out}")
            else:
                log(f"Failed to create subset file for cohort {cohort_name}", "ERROR")
        
        # 9. Quality check
        log(f"Quality check for {cohort_name} Multihetsep file:")
        check_cmd = f"head -n 10 {multihetsep_out}"
        run(check_cmd)
        
        log(f"Finished eSMC2 preparation for cohort: {cohort_name}")
    
    # Print summary
    log("eSMC2 preparation completed for all cohorts.")
    log("Important reminder:")
    log("  The Multihetsep file should not have a large number of '1's in the third column.")
    log("  Use 'less <multihetsep_file>' to inspect the files further.")
    
    # Print timing information
    end_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    log(f"Pipeline started at: {start_time}")
    log(f"Pipeline completed at: {end_time}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())