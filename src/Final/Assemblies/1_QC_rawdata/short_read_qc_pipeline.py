#!/usr/bin/env python3
"""
Short Read QC Pipeline

This script performs quality control processing for short read sequencing data, including:
- FastQC report generation for initial quality assessment
- Trimmomatic adapter trimming for reads with adapter contamination
- fastp-based polyG trimming and artifact removal
- Post-processing quality checks

The script processes multiple samples simultaneously and generates comprehensive reports
to track quality metrics throughout the pipeline.

Author: Max Borgmann
Date: March 2025
"""

import os
import subprocess
import logging
import argparse
import sys
import re
import json
import shutil
from datetime import datetime
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Set up logging
def setup_logging(log_directory, log_filename="short_read_qc.log"):
    """Configure logging to write to both file and console with detailed formatting."""
    os.makedirs(log_directory, exist_ok=True)
    log_file_path = os.path.join(log_directory, log_filename)
    
    # Create logger
    logger = logging.getLogger('qc_pipeline')
    logger.setLevel(logging.DEBUG)
    
    # Create file handler
    file_handler = logging.FileHandler(log_file_path)
    file_handler.setLevel(logging.DEBUG)
    
    # Create console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    
    # Create formatters
    file_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    console_formatter = logging.Formatter('%(message)s')
    
    # Add formatters to handlers
    file_handler.setFormatter(file_formatter)
    console_handler.setFormatter(console_formatter)
    
    # Add handlers to logger
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)
    
    return logger, log_file_path

# Utility function to run commands
def run(cmd, workdir=None, desc=None):
    """Execute a shell command and log the output."""
    if desc:
        logger.info(f"🔧 {desc}")
    logger.info(f"Running: {cmd}")
    
    try:
        process = subprocess.run(
            cmd, 
            shell=True, 
            check=True, 
            cwd=workdir,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True
        )
        logger.debug(f"Command stdout: {process.stdout.strip() if process.stdout else 'None'}")
        if process.stderr:
            logger.debug(f"Command stderr: {process.stderr.strip()}")
        return True, process.stdout
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed with exit code {e.returncode}")
        logger.error(f"Command stdout: {e.stdout.strip() if e.stdout else 'None'}")
        logger.error(f"Command stderr: {e.stderr.strip() if e.stderr else 'None'}")
        return False, e.stderr

# Function to check for adapter content in FastQC report
def check_adapter_content(fastqc_data_path):
    """Parse FastQC data file to check for adapter content."""
    adapter_content = False
    
    try:
        with open(fastqc_data_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line == ">>Adapter Content\tfail" or line == ">>Adapter Content\twarn":
                    adapter_content = True
                    break
    except Exception as e:
        logger.error(f"Error checking adapter content: {e}")
    
    return adapter_content

# Function to check for polyG content
def check_polyg_content(fastqc_data_path, threshold=10):
    """Parse FastQC data file to check for polyG content."""
    polyg_content = False
    in_over_represented = False
    
    try:
        with open(fastqc_data_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line == ">>Overrepresented sequences\tfail" or line == ">>Overrepresented sequences\twarn":
                    in_over_represented = True
                elif line == ">>END_MODULE":
                    in_over_represented = False
                elif in_over_represented and "GGGGGGGGGG" in line:
                    polyg_content = True
                    break
    except Exception as e:
        logger.error(f"Error checking polyG content: {e}")
    
    return polyg_content

# Function to extract basic sequence stats from FastQC data
def extract_fastqc_stats(fastqc_data_path):
    """Extract basic sequence statistics from FastQC data file."""
    stats = {
        "total_sequences": 0,
        "sequence_length": "0",
        "gc_content": 0
    }
    
    try:
        with open(fastqc_data_path, 'r') as f:
            in_basic_stats = False
            for line in f:
                line = line.strip()
                if line == ">>Basic Statistics":
                    in_basic_stats = True
                elif line == ">>END_MODULE":
                    in_basic_stats = False
                elif in_basic_stats:
                    if "Total Sequences" in line:
                        stats["total_sequences"] = int(line.split("\t")[1])
                    elif "Sequence length" in line:
                        stats["sequence_length"] = line.split("\t")[1]
                    elif "%GC" in line:
                        stats["gc_content"] = float(line.split("\t")[1])
    except Exception as e:
        logger.error(f"Error extracting FastQC stats: {e}")
    
    return stats

# Function to run FastQC
def run_fastqc(fastq_files, output_dir, threads):
    """Run FastQC on input FASTQ files."""
    os.makedirs(output_dir, exist_ok=True)
    
    # Group FASTQ files into a single command for efficiency
    fastq_list = " ".join(fastq_files)
    fastqc_cmd = f"fastqc -t {threads} -o {output_dir} {fastq_list}"
    
    success, _ = run(fastqc_cmd, desc="Running FastQC")
    
    # Parse FastQC results
    results = {}
    for fastq_file in fastq_files:
        basename = os.path.basename(fastq_file).replace(".fastq", "").replace(".fq", "")
        fastqc_data_path = os.path.join(output_dir, f"{basename}_fastqc", "fastqc_data.txt")
        
        if os.path.exists(fastqc_data_path):
            has_adapter = check_adapter_content(fastqc_data_path)
            has_polyg = check_polyg_content(fastqc_data_path)
            basic_stats = extract_fastqc_stats(fastqc_data_path)
            
            results[fastq_file] = {
                "has_adapter": has_adapter,
                "has_polyg": has_polyg,
                "fastqc_data_path": fastqc_data_path,
                "stats": basic_stats
            }
            
            logger.info(f"FastQC results for {basename}:")
            logger.info(f"  - Total sequences: {basic_stats['total_sequences']:,}")
            logger.info(f"  - Sequence length: {basic_stats['sequence_length']}")
            logger.info(f"  - GC content: {basic_stats['gc_content']}%")
            logger.info(f"  - Adapter content: {'Present' if has_adapter else 'Not detected'}")
            logger.info(f"  - PolyG content: {'Present' if has_polyg else 'Not detected'}")
        else:
            logger.warning(f"FastQC data not found for {basename}")
            results[fastq_file] = {
                "has_adapter": None,
                "has_polyg": None,
                "fastqc_data_path": None,
                "stats": {"total_sequences": 0, "sequence_length": "0", "gc_content": 0}
            }
    
    return results

# Function to run Trimmomatic
def run_trimmomatic(r1_file, r2_file, output_dir, adapter_file, threads):
    """Run Trimmomatic on paired FASTQ files."""
    os.makedirs(output_dir, exist_ok=True)
    
    # Define output files
    basename_r1 = os.path.basename(r1_file).replace(".fastq", "").replace(".fq", "")
    basename_r2 = os.path.basename(r2_file).replace(".fastq", "").replace(".fq", "")
    
    paired_r1 = os.path.join(output_dir, f"{basename_r1}.paired.fastq")
    unpaired_r1 = os.path.join(output_dir, f"{basename_r1}.unpaired.fastq")
    paired_r2 = os.path.join(output_dir, f"{basename_r2}.paired.fastq")
    unpaired_r2 = os.path.join(output_dir, f"{basename_r2}.unpaired.fastq")
    
    # Trimmomatic command
    trim_cmd = (
        f"trimmomatic PE -threads {threads} "
        f"{r1_file} {r2_file} "
        f"{paired_r1} {unpaired_r1} {paired_r2} {unpaired_r2} "
        f"ILLUMINACLIP:{adapter_file}:2:30:10 "
        f"LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
    )
    
    success, output = run(trim_cmd, desc=f"Running Trimmomatic on {basename_r1}/{basename_r2}")
    
    # Parse Trimmomatic output to get statistics
    stats = {
        "input_read_pairs": 0,
        "surviving_read_pairs": 0,
        "forward_only_surviving": 0,
        "reverse_only_surviving": 0,
        "dropped_read_pairs": 0
    }
    
    if success and output:
        # Try to parse Trimmomatic stats from output
        match = re.search(r'Input Read Pairs: (\d+).+Both Surviving: (\d+).+Forward Only Surviving: (\d+).+Reverse Only Surviving: (\d+).+Dropped: (\d+)', 
                          output, re.DOTALL)
        if match:
            stats["input_read_pairs"] = int(match.group(1))
            stats["surviving_read_pairs"] = int(match.group(2))
            stats["forward_only_surviving"] = int(match.group(3))
            stats["reverse_only_surviving"] = int(match.group(4))
            stats["dropped_read_pairs"] = int(match.group(5))
    
    logger.info(f"Trimmomatic processing results for {basename_r1}/{basename_r2}:")
    logger.info(f"  - Input read pairs: {stats['input_read_pairs']:,}")
    logger.info(f"  - Surviving read pairs: {stats['surviving_read_pairs']:,} ({stats['surviving_read_pairs']/stats['input_read_pairs']*100:.2f}% if stats['input_read_pairs'] > 0 else 0}}%)")
    logger.info(f"  - Forward only surviving: {stats['forward_only_surviving']:,}")
    logger.info(f"  - Reverse only surviving: {stats['reverse_only_surviving']:,}")
    logger.info(f"  - Dropped read pairs: {stats['dropped_read_pairs']:,}")
    
    return {
        "paired_r1": paired_r1,
        "unpaired_r1": unpaired_r1,
        "paired_r2": paired_r2,
        "unpaired_r2": unpaired_r2,
        "stats": stats
    }

# Function to run fastp
def run_fastp(r1_file, r2_file, output_dir, threads):
    """Run fastp for polyG trimming and artifact removal."""
    os.makedirs(output_dir, exist_ok=True)
    
    # Define output files
    basename_r1 = os.path.basename(r1_file).replace(".fastq", "").replace(".fq", "").replace(".paired", "")
    basename_r2 = os.path.basename(r2_file).replace(".fastq", "").replace(".fq", "").replace(".paired", "")
    
    out_r1 = os.path.join(output_dir, f"{basename_r1}.fastp.fastq")
    out_r2 = os.path.join(output_dir, f"{basename_r2}.fastp.fastq")
    json_report = os.path.join(output_dir, f"{basename_r1}_{basename_r2}.fastp.json")
    html_report = os.path.join(output_dir, f"{basename_r1}_{basename_r2}.fastp.html")
    
    # fastp command with polyG trimming
    fastp_cmd = (
        f"fastp --in1 {r1_file} --in2 {r2_file} "
        f"--out1 {out_r1} --out2 {out_r2} "
        f"--json {json_report} --html {html_report} "
        f"--thread {threads} --trim_poly_g --trim_poly_x "
        f"--detect_adapter_for_pe --cut_right --cut_right_window_size 4 "
        f"--cut_right_mean_quality 20 --correction"
    )
    
    success, _ = run(fastp_cmd, desc=f"Running fastp on {basename_r1}/{basename_r2}")
    
    # Parse fastp JSON report to get statistics
    stats = {
        "before_total_reads": 0,
        "before_total_bases": 0,
        "after_total_reads": 0,
        "after_total_bases": 0,
        "polyg_trimmed_reads": 0,
        "retained_reads_percent": 0,
        "retained_bases_percent": 0
    }
    
    if os.path.exists(json_report):
        try:
            with open(json_report, 'r') as f:
                report_data = json.load(f)
                before = report_data.get('summary', {}).get('before_filtering', {})
                after = report_data.get('summary', {}).get('after_filtering', {})
                
                stats["before_total_reads"] = before.get('total_reads', 0)
                stats["before_total_bases"] = before.get('total_bases', 0)
                stats["after_total_reads"] = after.get('total_reads', 0)
                stats["after_total_bases"] = after.get('total_bases', 0)
                stats["polyg_trimmed_reads"] = report_data.get('filtering_result', {}).get('poly_g_trimmed_reads', 0)
                
                if stats["before_total_reads"] > 0:
                    stats["retained_reads_percent"] = (stats["after_total_reads"] / stats["before_total_reads"]) * 100
                if stats["before_total_bases"] > 0:
                    stats["retained_bases_percent"] = (stats["after_total_bases"] / stats["before_total_bases"]) * 100
        except Exception as e:
            logger.error(f"Error parsing fastp JSON report: {e}")
    
    logger.info(f"fastp processing results for {basename_r1}/{basename_r2}:")
    logger.info(f"  - Input reads: {stats['before_total_reads']:,}")
    logger.info(f"  - Output reads: {stats['after_total_reads']:,}")
    logger.info(f"  - Retained reads: {stats['retained_reads_percent']:.2f}%")
    logger.info(f"  - Input bases: {stats['before_total_bases']:,}")
    logger.info(f"  - Output bases: {stats['after_total_bases']:,}")
    logger.info(f"  - Retained bases: {stats['retained_bases_percent']:.2f}%")
    logger.info(f"  - PolyG trimmed reads: {stats['polyg_trimmed_reads']:,}")
    
    return {
        "out_r1": out_r1,
        "out_r2": out_r2,
        "json_report": json_report,
        "html_report": html_report,
        "stats": stats
    }

# Function to create QC summary
def create_qc_summary(samples, initial_fastqc, trimmomatic_results, fastp_results, summary_path):
    """Create a summary of QC processing results."""
    summary = ["# Short Read QC Pipeline Summary\n"]
    summary.append(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    
    summary.append("## Samples Processed\n")
    for sample, files in samples.items():
        summary.append(f"- **{sample}**")
        summary.append(f"  - R1: {os.path.basename(files['r1'])}")
        summary.append(f"  - R2: {os.path.basename(files['r2'])}")
    summary.append("\n")
    
    summary.append("## Initial FastQC Results\n")
    summary.append("| Sample | File | Total Sequences | Sequence Length | GC% | Adapter Content | PolyG Content |")
    summary.append("|--------|------|----------------|----------------|-----|----------------|--------------|")
    for sample, files in samples.items():
        r1_result = initial_fastqc.get(files['r1'], {})
        r2_result = initial_fastqc.get(files['r2'], {})
        
        r1_stats = r1_result.get('stats', {"total_sequences": 0, "sequence_length": "0", "gc_content": 0})
        r2_stats = r2_result.get('stats', {"total_sequences": 0, "sequence_length": "0", "gc_content": 0})
        
        summary.append(f"| {sample} | R1 | {r1_stats['total_sequences']:,} | {r1_stats['sequence_length']} | {r1_stats['gc_content']}% | {'Yes' if r1_result.get('has_adapter') else 'No'} | {'Yes' if r1_result.get('has_polyg') else 'No'} |")
        summary.append(f"| {sample} | R2 | {r2_stats['total_sequences']:,} | {r2_stats['sequence_length']} | {r2_stats['gc_content']}% | {'Yes' if r2_result.get('has_adapter') else 'No'} | {'Yes' if r2_result.get('has_polyg') else 'No'} |")
    summary.append("\n")
    
    if trimmomatic_results:
        summary.append("## Trimmomatic Results\n")
        summary.append("| Sample | Input Pairs | Surviving Pairs | % Retained | Forward Only | Reverse Only | Dropped |")
        summary.append("|--------|------------|----------------|------------|--------------|--------------|---------|")
        for sample, result in trimmomatic_results.items():
            stats = result.get('stats', {})
            input_pairs = stats.get('input_read_pairs', 0)
            surv_pairs = stats.get('surviving_read_pairs', 0)
            percent = (surv_pairs / input_pairs * 100) if input_pairs > 0 else 0
            
            summary.append(f"| {sample} | {input_pairs:,} | {surv_pairs:,} | {percent:.2f}% | "
                          f"{stats.get('forward_only_surviving', 0):,} | {stats.get('reverse_only_surviving', 0):,} | "
                          f"{stats.get('dropped_read_pairs', 0):,} |")
        summary.append("\n")
    
    if fastp_results:
        summary.append("## fastp Processing Results\n")
        summary.append("| Sample | Input Reads | Output Reads | % Retained | Input Bases | Output Bases | % Retained | PolyG Trimmed |")
        summary.append("|--------|------------|--------------|------------|-------------|--------------|------------|---------------|")
        for sample, result in fastp_results.items():
            stats = result.get('stats', {})
            in_reads = stats.get('before_total_reads', 0)
            out_reads = stats.get('after_total_reads', 0)
            reads_pct = stats.get('retained_reads_percent', 0)
            
            in_bases = stats.get('before_total_bases', 0)
            out_bases = stats.get('after_total_bases', 0)
            bases_pct = stats.get('retained_bases_percent', 0)
            
            summary.append(f"| {sample} | {in_reads:,} | {out_reads:,} | {reads_pct:.2f}% | "
                          f"{in_bases:,} | {out_bases:,} | {bases_pct:.2f}% | "
                          f"{stats.get('polyg_trimmed_reads', 0):,} |")
        summary.append("\n")
    
    summary.append("## Pipeline Steps Performed\n")
    summary.append("1. Initial FastQC analysis")
    if trimmomatic_results:
        summary.append("2. Trimmomatic adapter trimming")
    summary.append(f"{'3' if trimmomatic_results else '2'}. fastp polyG trimming and quality filtering")
    summary.append(f"{'4' if trimmomatic_results else '3'}. Final quality assessment")
    
    # Write summary to file
    with open(summary_path, 'w') as f:
        f.write("\n".join(summary))
    
    logger.info(f"QC summary saved to: {summary_path}")
    return summary

# Function to generate quality plots
def generate_quality_plots(samples, fastqc_dir, trimmomatic_dir, fastp_dir, final_fastqc_dir, plot_dir):
    """Generate summary plots comparing quality metrics before and after processing."""
    os.makedirs(plot_dir, exist_ok=True)
    
    # Plot total reads before and after processing
    try:
        sample_names = []
        initial_reads = []
        final_reads = []
        
        for sample, _ in samples.items():
            # Get initial and final read counts for each sample
            # This would need to be implemented by parsing FastQC results
            pass
        
        # Create the plot
        # This would be implemented using matplotlib
        pass
        
    except Exception as e:
        logger.error(f"Error generating read count plot: {e}")
    
    # Additional plots could be implemented here
    
    return True

# Main function
def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Short Read QC Pipeline")
    
    parser.add_argument("-i", "--input_dir", 
                        help="Input directory containing FASTQ files")
    parser.add_argument("-f", "--input_files", nargs="+",
                        help="Input FASTQ files (specify both R1 and R2 for each sample)")
    parser.add_argument("-o", "--output_dir", default="./QC_results",
                        help="Output directory for QC results")
    parser.add_argument("-a", "--adapter_file", default="/usr/share/trimmomatic/adapters/TruSeq3-PE.fa",
                        help="Adapter file for Trimmomatic")
    parser.add_argument("-t", "--threads", type=int, default=8,
                        help="Number of CPU threads to use")
    parser.add_argument("--skip_trimmomatic", action="store_true",
                        help="Skip Trimmomatic even if adapter content is detected")
    parser.add_argument("--skip_fastp", action="store_true",
                        help="Skip fastp polyG trimming")
    parser.add_argument("--overwrite", action="store_true",
                        help="Overwrite existing output files")
    
    args = parser.parse_args()
    
    # Check if either input_dir or input_files is provided
    if not args.input_dir and not args.input_files:
        print("Error: Either --input_dir or --input_files must be provided")
        parser.print_help()
        sys.exit(1)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Set up logging
    global logger
    logger, log_file_path = setup_logging(args.output_dir)
    
    logger.info("=" * 60)
    logger.info("Short Read QC Pipeline")
    logger.info("=" * 60)
    logger.info(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Output directory: {os.path.abspath(args.output_dir)}")
    logger.info(f"Log file: {log_file_path}")
    logger.info(f"Threads: {args.threads}")
    
    # Create subdirectories
    fastqc_dir = os.path.join(args.output_dir, "01_fastqc")
    trimmomatic_dir = os.path.join(args.output_dir, "02_trimmomatic")
    fastp_dir = os.path.join(args.output_dir, "03_fastp")
    final_dir = os.path.join(args.output_dir, "04_cleaned_reads")
    final_fastqc_dir = os.path.join(args.output_dir, "05_final_fastqc")
    plot_dir = os.path.join(args.output_dir, "06_quality_plots")
    
    os.makedirs(fastqc_dir, exist_ok=True)
    os.makedirs(trimmomatic_dir, exist_ok=True)
    os.makedirs(fastp_dir, exist_ok=True)
    os.makedirs(final_dir, exist_ok=True)
    os.makedirs(final_fastqc_dir, exist_ok=True)
    os.makedirs(plot_dir, exist_ok=True)
    
    # Collect input files
    input_files = []
    if args.input_files:
        input_files = args.input_files
    elif args.input_dir:
        for file in os.listdir(args.input_dir):
            if file.endswith(".fastq") or file.endswith(".fq"):
                input_files.append(os.path.join(args.input_dir, file))
    
    logger.info(f"Found {len(input_files)} input files")
    
    # Group input files by sample
    samples = {}
    r1_files = [f for f in input_files if "_R1" in f]
    r2_files = [f for f in input_files if "_R2" in f]
    
    if len(r1_files) != len(r2_files):
        logger.error("Number of R1 files does not match number of R2 files")
        sys.exit(1)
    
    for r1 in r1_files:
        sample_name = os.path.basename(r1).split("_R1")[0]
        r2 = None
        
        # Find matching R2 file
        for potential_r2 in r2_files:
            if os.path.basename(potential_r2).split("_R2")[0] == sample_name:
                r2 = potential_r2
                break
        
        if not r2:
            logger.error(f"Cannot find matching R2 file for {r1}")
            sys.exit(1)
        
        samples[sample_name] = {
            "r1": r1,
            "r2": r2
        }
    
    logger.info(f"Found {len(samples)} sample(s):")
    for sample, files in samples.items():
        logger.info(f"  - {sample}: {os.path.basename(files['r1'])} & {os.path.basename(files['r2'])}")
    
    # Step 1: Run initial FastQC
    logger.info("\n" + "=" * 60)
    logger.info("Step 1: Running initial FastQC analysis")
    logger.info("=" * 60)
    
    all_fastq_files = [files['r1'] for _, files in samples.items()] + [files['r2'] for _, files in samples.items()]
    initial_fastqc = run_fastqc(all_fastq_files, fastqc_dir, args.threads)
    
    # Check if any files need adapter trimming
    need_adapter_trimming = any(result.get('has_adapter', False) for result in initial_fastqc.values())
    
    # Step 2: Run Trimmomatic if needed and not skipped
    trimmomatic_results = {}
    if need_adapter_trimming and not args.skip_trimmomatic:
        logger.info("\n" + "=" * 60)
        logger.info("Step 2: Running Trimmomatic for adapter removal")
        logger.info("=" * 60)
        
        for sample, files in samples.items():
            r1_result = initial_fastqc.get(files['r1'], {})
            r2_result = initial_fastqc.get(files['r2'], {})
            
            if r1_result.get('has_adapter', False) or r2_result.get('has_adapter', False):
                logger.info(f"Running Trimmomatic on {sample} (adapter content detected)")
                result = run_trimmomatic(files['r1'], files['r2'], trimmomatic_dir, args.adapter_file, args.threads)
                trimmomatic_results[sample] = result
                
                # Update files to use trimmed files for fastp
                samples[sample]['r1'] = result['paired_r1']
                samples[sample]['r2'] = result['paired_r2']
            else:
                logger.info(f"Skipping Trimmomatic for {sample} (no adapter content detected)")
    elif args.skip_trimmomatic:
        logger.info("Trimmomatic step skipped by user request")
    else:
        logger.info("No adapter content detected in any files, skipping Trimmomatic")
    
    # Step 3: Run fastp for polyG trimming if not skipped
    fastp_results = {}
    if not args.skip_fastp:
        logger.info("\n" + "=" * 60)
        logger.info("Step 3: Running fastp for polyG trimming and quality filtering")
        logger.info("=" * 60)
        
        for sample, files in samples.items():
            logger.info(f"Running fastp on {sample}")
            result = run_fastp(files['r1'], files['r2'], fastp_dir, args.threads)
            fastp_results[sample] = result
            
            # Copy final cleaned files to the final directory
            final_r1 = os.path.join(final_dir, f"{sample}_R1.clean.fastq")
            final_r2 = os.path.join(final_dir, f"{sample}_R2.clean.fastq")
            
            shutil.copy2(result['out_r1'], final_r1)
            shutil.copy2(result['out_r2'], final_r2)
            
            logger.info(f"Final cleaned files copied to:")
            logger.info(f"  - {final_r1}")
            logger.info(f"  - {final_r2}")
    else:
        logger.info("fastp step skipped by user request")
        
        # If fastp is skipped, copy the current files (original or trimmomatic-processed) to final directory
        for sample, files in samples.items():
            final_r1 = os.path.join(final_dir, f"{sample}_R1.clean.fastq")
            final_r2 = os.path.join(final_dir, f"{sample}_R2.clean.fastq")
            
            shutil.copy2(files['r1'], final_r1)
            shutil.copy2(files['r2'], final_r2)
            
            logger.info(f"Files copied to final directory for {sample}:")
            logger.info(f"  - {final_r1}")
            logger.info(f"  - {final_r2}")
    
    # Step 4: Run final FastQC on cleaned reads
    logger.info("\n" + "=" * 60)
    logger.info("Step 4: Running final FastQC on cleaned reads")
    logger.info("=" * 60)
    
    final_fastq_files = [os.path.join(final_dir, f"{sample}_R1.clean.fastq") for sample in samples.keys()]
    final_fastq_files += [os.path.join(final_dir, f"{sample}_R2.clean.fastq") for sample in samples.keys()]
    
    final_fastqc = run_fastqc(final_fastq_files, final_fastqc_dir, args.threads)
    
    # Step 5: Generate quality plots
    logger.info("\n" + "=" * 60)
    logger.info("Step 5: Generating quality comparison plots")
    logger.info("=" * 60)
    
    generate_quality_plots(samples, fastqc_dir, trimmomatic_dir, fastp_dir, final_fastqc_dir, plot_dir)
    
    # Create summary report
    logger.info("\n" + "=" * 60)
    logger.info("Creating QC summary report")
    logger.info("=" * 60)
    
    summary_path = os.path.join(args.output_dir, "QC_summary.md")
    summary = create_qc_summary(samples, initial_fastqc, trimmomatic_results, fastp_results, summary_path)
    
    logger.info("\n" + "=" * 60)
    logger.info("Short Read QC Pipeline Completed Successfully")
    logger.info("=" * 60)
    logger.info(f"End time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"QC summary: {summary_path}")
    logger.info(f"Cleaned reads directory: {final_dir}")

if __name__ == "__main__":
    main()