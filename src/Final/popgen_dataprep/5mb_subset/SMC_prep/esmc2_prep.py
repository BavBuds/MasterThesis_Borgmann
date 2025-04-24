import os
import subprocess
import sys
import argparse
import shutil
import glob
import re
import time
import json
from datetime import datetime

def check_dependencies():
    """Check if required software is installed."""
    required_tools = {
        "gatk": "GATK",
        "bcftools": "Bcftools",
        "vcf2bed": "Bedops",
        "bedtools": "Bedtools",
        "samtools": "Samtools",
        "wget": "Wget",
        "bwa-mem2": "BWA-MEM2"
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
    parser.add_argument("--depth-filter", type=int, default=10, 
                        help="Minimum depth filter (default: 10)")
    parser.add_argument("--memory", default="16g", 
                        help="Memory allocation for Java (default: 16g)")
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
    parser.add_argument("--use-existing-cohorts", action="store_true",
                        help="Use existing cohorts found in output directory instead of creating new ones")
    parser.add_argument("--skip-chromosome-splitting", action="store_true",
                        help="Skip processing by chromosome and process whole VCF at once")
 

    
    # Add remapping arguments
    parser.add_argument("--remap-sample", 
                        help="Sample to remap (e.g., 'HT2')")
    parser.add_argument("--remap-fastq1", 
                        help="Path to first FASTQ file for remapping")
    parser.add_argument("--remap-fastq2", 
                        help="Path to second FASTQ file for remapping")
    parser.add_argument("--remap-to-reference", 
                        help="Path to reference FASTA for remapping")
    parser.add_argument("--remap-output-dir", 
                        help="Directory for remapped output (default: <output-dir>/<cohort-name>)")
    parser.add_argument('--override-ref',  action='append',
                    help='sample_id:/full/path/to/new_reference.fasta')
    parser.add_argument('--override-bam',  action='append',
                    help='sample_id:/full/path/to/filtered_or_fixed.bam')

    
    return parser.parse_args()


def parse_override_args(override_list):
    """
    Parse override arguments in the format 'sample_id:/path/to/file'.
    
    Args:
        override_list (list): List of override arguments
        
    Returns:
        dict: Dictionary of {sample_id: file_path} overrides
    """
    overrides = {}
    if not override_list:
        return overrides
        
    for override in override_list:
        try:
            sample_id, file_path = override.split(':', 1)
            sample_id = sample_id.strip()
            file_path = file_path.strip()
            
            if not sample_id or not file_path:
                log(f"Invalid override format: {override}. Expected 'sample_id:/path/to/file'", "WARNING")
                continue
                
            if not os.path.exists(file_path):
                log(f"Warning: File does not exist: {file_path}", "WARNING")
                
            overrides[sample_id] = file_path
            log(f"Override set for {sample_id}: {file_path}", "INFO")
        except ValueError:
            log(f"Invalid override format: {override}. Expected 'sample_id:/path/to/file'", "WARNING")
            
    return overrides


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

def index_vcf(vcf_path, threads=1):
    """
    Index a VCF file using bcftools.
    
    Args:
        vcf_path (str): Path to the VCF file
        threads (int): Number of threads to use
        
    Returns:
        bool: True if indexing was successful, False otherwise
    """
    if not os.path.exists(vcf_path):
        log(f"Cannot index VCF: File not found - {vcf_path}", "ERROR")
        return False
    
    # Check if index already exists
    if os.path.exists(vcf_path + ".tbi") or os.path.exists(vcf_path + ".csi"):
        log(f"VCF index already exists for {vcf_path}", "INFO")
        return True
    
    # Index the VCF
    log(f"Indexing VCF: {vcf_path}", "INFO")
    index_cmd = f"bcftools index --threads {threads} {vcf_path}"
    return run(index_cmd)[0]

def patch_generate_multihetsep(script_path):
    log("Patching generate_multihetsep.py script to handle unsorted positions")

    try:
        with open(script_path, "r") as f:
            script = f.read()
    except Exception as e:
        log(f"Cannot read script: {e}", "ERROR")
        return False

    if "if pos < self.lastPos" in script:
        log("Patch already present – nothing to do.", "INFO")
        return True          # ← let the pipeline proceed

    if "assert pos >= self.lastPos" in script:
        patched = script.replace(
            "assert pos >= self.lastPos",
            "if pos < self.lastPos: return False  # Skip positions that go backwards"
        )
        try:
            with open(script_path, "w") as f:
                f.write(patched)
            log("Successfully patched generate_multihetsep.py", "INFO")
            return True
        except Exception as e:
            log(f"Write failed: {e}", "ERROR")
            return False

    log("Unknown script structure – skipping patch but continuing.", "INFO")
    return True


def remap_reads(fastq1, fastq2, reference_fasta, sample_name, output_dir, threads):
    """
    Remap reads to a reference genome using BWA-MEM2.
    
    Args:
        fastq1 (str): Path to first FASTQ file
        fastq2 (str): Path to second FASTQ file
        reference_fasta (str): Path to reference FASTA
        sample_name (str): Sample name
        output_dir (str): Output directory
        threads (int): Number of threads
    
    Returns:
        str: Path to remapped BAM file with read groups
    """
    log(f"Remapping {sample_name} reads to reference: {reference_fasta}", "INFO")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    fixed_bams_dir = os.path.join(output_dir, "fixed_bams")
    os.makedirs(fixed_bams_dir, exist_ok=True)
    
    # Output paths
    sorted_bam = os.path.join(output_dir, f"{sample_name}_sorted.bam")
    fixed_bam = os.path.join(fixed_bams_dir, f"{sample_name}_filtered.with_rg.bam")
    
    # Skip if output already exists
    if os.path.exists(fixed_bam):
        log(f"Remapped BAM already exists: {fixed_bam}", "INFO")
        return fixed_bam
    
    # Index reference if needed
    bwa_index = reference_fasta + ".bwt.2bit.64"
    if not os.path.exists(bwa_index):
        log(f"Indexing reference with BWA-MEM2: {reference_fasta}", "INFO")
        index_cmd = f"bwa-mem2 index {reference_fasta}"
        if not run(index_cmd)[0]:
            log(f"Failed to index reference: {reference_fasta}", "ERROR")
            return None
    
    # Align reads
    log(f"Aligning {sample_name} reads to reference using BWA-MEM2", "INFO")
    align_cmd = (
        f"bwa-mem2 mem -t {threads} {reference_fasta} {fastq1} {fastq2} | "
        f"samtools view -bS - | "
        f"samtools sort -@ {threads} -o {sorted_bam} -"
    )
    if not run(align_cmd)[0]:
        log(f"Failed to align {sample_name} reads to reference", "ERROR")
        return None
    
    # Add read groups
    log(f"Adding read groups to {sample_name} BAM", "INFO")
    rg_cmd = (
        f"samtools addreplacerg -r '@RG\\tID:{sample_name}\\tSM:{sample_name}\\tPL:ILLUMINA' "
        f"-o {fixed_bam} {sorted_bam}"
    )
    if not run(rg_cmd)[0]:
        log(f"Failed to add read groups to {sample_name} BAM", "ERROR")
        return sorted_bam
    
    # Index fixed BAM
    log(f"Indexing {sample_name} BAM", "INFO")
    index_cmd = f"samtools index {fixed_bam}"
    if not run(index_cmd)[0]:
        log(f"Failed to index {sample_name} BAM", "ERROR")
        return fixed_bam
    
    log(f"Successfully remapped {sample_name} reads to reference: {fixed_bam}", "INFO")
    return fixed_bam

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

def discover_existing_cohorts(output_dir):
    """
    Discover existing cohorts in the output directory by looking for cohort metadata files.
    
    Returns:
        dict: Dictionary of cohort information keyed by cohort name.
    """
    existing_cohorts = {}
    
    # Check for cohort directories
    for cohort_dir in glob.glob(os.path.join(output_dir, "*")):
        if not os.path.isdir(cohort_dir):
            continue
            
        cohort_name = os.path.basename(cohort_dir)
        if cohort_name == "GVCFs":  # Skip GVCFs directory
            continue
            
        metadata_file = os.path.join(cohort_dir, "cohort_metadata.json")
        if os.path.exists(metadata_file):
            # Found a cohort metadata file, load it
            try:
                with open(metadata_file, 'r') as f:
                    cohort_info = json.load(f)
                    log(f"Found existing cohort: {cohort_name}", "INFO")
                    existing_cohorts[cohort_name] = cohort_info
            except json.JSONDecodeError:
                log(f"Failed to parse cohort metadata file: {metadata_file}", "WARNING")
        else:
            # Check if this is a cohort directory by looking for common cohort files
            cohort_files = [
                f"{cohort_name}.5mb_subset.g.vcf.gz",
                f"{cohort_name}.wholeGenome.g.vcf.gz",
                f"{cohort_name}.allsites.geno.vcf.gz"
            ]
            
            for cf in cohort_files:
                if os.path.exists(os.path.join(cohort_dir, cf)):
                    log(f"Found existing cohort without metadata: {cohort_name}", "INFO")
                    
                    # Try to infer cohort information by listing sample VCFs
                    sample_vcfs = glob.glob(os.path.join(cohort_dir, f"{cohort_name}.*.vcf.gz"))
                    sample_names = []
                    for vcf in sample_vcfs:
                        match = re.search(rf"{cohort_name}\.([^\.]+)\.vcf\.gz$", vcf)
                        if match and match.group(1) not in ["final", "allsites"]:
                            sample_names.append(match.group(1))
                    
                    if sample_names:
                        log(f"Inferred samples for cohort {cohort_name}: {', '.join(sample_names)}", "INFO")
                        
                        # Create a best-guess metadata file
                        cohort_info = {
                            "name": cohort_name,
                            "samples": sample_names,
                            # Using HT reference by default for existing cohorts
                            "reference": "/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/Popgen_analysis/dataprep/5mb_subset/Input_fastas/HT/HT_assembly_5mb_subset.fasta",
                            "created": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                            "inferred": True
                        }
                        
                        try:
                            with open(metadata_file, 'w') as f:
                                json.dump(cohort_info, f, indent=2)
                            log(f"Created metadata file for existing cohort: {cohort_name}", "INFO")
                            existing_cohorts[cohort_name] = cohort_info
                        except Exception as e:
                            log(f"Failed to create metadata file for cohort {cohort_name}: {str(e)}", "WARNING")
                    
                    break
    
    return existing_cohorts

def save_cohort_metadata(cohort_dir, cohort_name, sample_list, reference_fasta):
    """
    Save cohort metadata to a JSON file.
    
    Args:
        cohort_dir (str): Directory where the cohort is processed
        cohort_name (str): Name of the cohort
        sample_list (list): List of (sample_name, bam_path) tuples
        reference_fasta (str): Path to the reference FASTA used for this cohort
    """
    metadata_file = os.path.join(cohort_dir, "cohort_metadata.json")
    
    # Create metadata dictionary
    metadata = {
        "name": cohort_name,
        "samples": [s for s, _ in sample_list],
        "reference": reference_fasta,
        "created": datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    }
    
    # Write to file
    try:
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)
        log(f"Saved cohort metadata to {metadata_file}", "INFO")
        return True
    except Exception as e:
        log(f"Failed to save cohort metadata: {str(e)}", "ERROR")
        return False

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

def find_gvcf_files(gvcf_dir, region_name="*"):
    """
    Find GVCF files in the GVCF directory.
    
    Args:
        gvcf_dir (str): Directory containing GVCF files
        region_name (str): Region name pattern to match (default: "*" to match any)
        
    Returns:
        dict: Dictionary of GVCF files keyed by sample name
    """
    gvcfs = {}
    
    # Look for sample.region.g.vcf.gz files
    for gvcf_file in glob.glob(os.path.join(gvcf_dir, f"*.{region_name}.g.vcf.gz")):
        # Extract sample name from filename
        basename = os.path.basename(gvcf_file)
        sample_name = basename.split('.')[0]  # Assuming format is sample.region.g.vcf.gz
        
        gvcfs[sample_name] = gvcf_file
        log(f"Found GVCF for {sample_name}: {gvcf_file}", "INFO")
    
    return gvcfs

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
    log("Configuration:")
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
    if args.remap_sample:
        log(f"  Remapping sample: {args.remap_sample}")
        log(f"  Remapping FASTQ1: {args.remap_fastq1}")
        log(f"  Remapping FASTQ2: {args.remap_fastq2}")
        log(f"  Remapping to reference: {args.remap_to_reference}")
    if args.use_existing_cohorts:
        log("  Using existing cohorts if found")
    if args.skip_chromosome_splitting:
        log("  Skipping chromosome-by-chromosome processing")
    if args.override_ref:
        log(f"  Reference overrides requested: {len(args.override_ref)}")
    if args.override_bam:
        log(f"  BAM overrides requested: {len(args.override_bam)}")
    
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
    
    # Apply reference fasta overrides
    if args.override_ref:
        ref_overrides = parse_override_args(args.override_ref)
        for sample_id, ref_path in ref_overrides.items():
            if os.path.exists(ref_path):
                reference_fastas[sample_id] = ref_path
                log(f"Overriding reference for {sample_id} with {ref_path}", "INFO")
                
                # Make sure the override reference is indexed
                if index_reference_fasta(ref_path):
                    log(f"Successfully indexed override reference for {sample_id}", "INFO")
                else:
                    log(f"Failed to index override reference for {sample_id}", "WARNING")
            else:
                log(f"Override reference file does not exist: {ref_path}", "ERROR")
    
    # Apply BAM file overrides
    if args.override_bam:
        bam_overrides = parse_override_args(args.override_bam)
        for sample_id, bam_path in bam_overrides.items():
            if os.path.exists(bam_path):
                filtered_bams[sample_id] = bam_path
                log(f"Overriding BAM for {sample_id} with {bam_path}", "INFO")
            else:
                log(f"Override BAM file does not exist: {bam_path}", "ERROR")
    
    # Handle remapping if requested
    if args.remap_sample and args.remap_fastq1 and args.remap_fastq2 and args.remap_to_reference:
        log(f"Remapping {args.remap_sample} to {args.remap_to_reference}")
        
        # Set default output directory based on reference name
        if not args.remap_output_dir:
            # Extract cohort name from reference path (assuming it contains the sample name)
            ref_basename = os.path.basename(os.path.dirname(args.remap_to_reference))
            args.remap_output_dir = os.path.join(args.output_dir, f"{ref_basename}tot")
        
        # Run remapping
        remapped_bam = remap_reads(
            args.remap_fastq1,
            args.remap_fastq2,
            args.remap_to_reference,
            args.remap_sample,
            args.remap_output_dir,
            args.threads
        )
        
        if remapped_bam:
            log(f"Successfully remapped {args.remap_sample} to {args.remap_to_reference}")
            
            # Update the filtered_bams dictionary
            filtered_bams[args.remap_sample] = remapped_bam
            log(f"Using remapped BAM for {args.remap_sample} in further processing")
        else:
            log(f"Failed to remap {args.remap_sample}", "ERROR")
    
    # Check and fix read groups in BAM files
    fixed_bams = {}
    sample_names = {}
    for sample_name, bam_file in filtered_bams.items():
        fixed_bam, actual_sample = check_and_fix_read_groups(bam_file, sample_name)
        fixed_bams[sample_name] = fixed_bam
        sample_names[sample_name] = actual_sample
    
    # Discover existing cohorts
    existing_cohorts = {}
    if args.use_existing_cohorts:
        existing_cohorts = discover_existing_cohorts(args.output_dir)
        log(f"Found {len(existing_cohorts)} existing cohorts")
    
    # Prompt user to define cohorts if no existing cohorts or not using existing cohorts
    cohorts = existing_cohorts.copy()
    if not args.use_existing_cohorts or not existing_cohorts:
        user_cohorts = prompt_user_for_cohorts(filtered_bams, reference_fastas)
        cohorts.update(user_cohorts)
    
    if not cohorts:
        log("No cohorts defined. Exiting.", "WARNING")
        return 1
    
    # Region parameter for GATK
    region_param = f"-L {args.region}" if args.region else ""
    region_name = args.region.replace(":", "_") if args.region else "5mb_subset"
    
    # Find existing GVCFs - look for both "5mb_subset" and "wholeGenome" names
    gvcfs = find_gvcf_files(gvcf_dir, region_name)
    wholegenome_gvcfs = find_gvcf_files(gvcf_dir, "wholeGenome")
    
    # Merge both dictionaries, preferring the region-specific GVCFs
    all_gvcfs = {**wholegenome_gvcfs, **gvcfs}
    
    # Step 1: Run HaplotypeCaller on each BAM file
    for sample_name, bam_file in fixed_bams.items():
        # Skip if GVCF already exists in either format
        if sample_name in all_gvcfs:
            log(f"GVCF already exists for {sample_name}: {all_gvcfs[sample_name]}")
            continue
            
        # Special case for remapped sample - use the reference specified for remapping
        if sample_name == args.remap_sample and args.remap_to_reference:
            reference_fasta = args.remap_to_reference
            log(f"Using specified remapping reference for {sample_name}: {reference_fasta}")
        elif sample_name not in reference_fastas:
            log(f"Skipping HaplotypeCaller for {sample_name} - no reference fasta found", "WARNING")
            continue
        else:
            reference_fasta = reference_fastas[sample_name]
            
        actual_sample = sample_names[sample_name]
        log(f"Running HaplotypeCaller on sample: {sample_name} (using sample name: {actual_sample})")
        
        # Check for both naming patterns
        gvcf_out = os.path.join(gvcf_dir, f"{sample_name}.{region_name}.g.vcf.gz")
        wholegenome_gvcf_out = os.path.join(gvcf_dir, f"{sample_name}.wholeGenome.g.vcf.gz")
        
        if os.path.exists(gvcf_out):
            log(f"GVCF already exists for {sample_name}: {gvcf_out}")
            all_gvcfs[sample_name] = gvcf_out
            continue
        elif os.path.exists(wholegenome_gvcf_out):
            log(f"GVCF already exists for {sample_name}: {wholegenome_gvcf_out}")
            all_gvcfs[sample_name] = wholegenome_gvcf_out
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
            # Index the newly created GVCF
            if index_vcf(gvcf_out, args.threads):
                all_gvcfs[sample_name] = gvcf_out
                log(f"Successfully created and indexed GVCF for {sample_name}")
            else:
                log(f"Created GVCF but failed to index it for {sample_name}", "WARNING")
                all_gvcfs[sample_name] = gvcf_out
        else:
            log(f"Failed to create GVCF for {sample_name}", "ERROR")
    
    if not all_gvcfs:
        log("No GVCFs were created or found. Exiting.", "ERROR")
        return 1
    
    # Process each cohort
    for cohort_name, cohort_info in cohorts.items():
        log(f"Processing cohort: {cohort_name}")
        cohort_dir = os.path.join(args.output_dir, cohort_name)
        os.makedirs(cohort_dir, exist_ok=True)
        
        # Get samples and reference for this cohort
        if "samples" in cohort_info and isinstance(cohort_info["samples"], list):
            if isinstance(cohort_info["samples"][0], tuple):
                # Format from prompt_user_for_cohorts: [(sample_name, bam_path), ...]
                sample_list = cohort_info["samples"]
            else:
                # Format from metadata file: [sample_name, ...]
                sample_list = [(s, filtered_bams.get(s, "")) for s in cohort_info["samples"]]
        else:
            log(f"Invalid sample list for cohort {cohort_name}", "ERROR")
            continue
            
        reference_fasta = cohort_info.get("reference")
        if not reference_fasta or not os.path.exists(reference_fasta):
            log(f"Invalid reference fasta for cohort {cohort_name}: {reference_fasta}", "ERROR")
            continue
        
        # Save metadata if it doesn't exist
        metadata_file = os.path.join(cohort_dir, "cohort_metadata.json")
        if not os.path.exists(metadata_file):
            save_cohort_metadata(cohort_dir, cohort_name, sample_list, reference_fasta)
        
        # Get GVCF files for samples in this cohort
        cohort_gvcfs = []
        missing_samples = []
        for sample_name, _ in sample_list:
            if sample_name in all_gvcfs:
                cohort_gvcfs.append((sample_name, all_gvcfs[sample_name]))
            else:
                missing_samples.append(sample_name)
        
        if missing_samples:
            log(f"Missing GVCFs for samples in cohort {cohort_name}: {', '.join(missing_samples)}", "WARNING")
            if not cohort_gvcfs:
                log(f"Skipping cohort {cohort_name} - no GVCFs available", "ERROR")
                continue
        
        # 1. Combine GVCFs
        combined_gvcf = os.path.join(cohort_dir, f"{cohort_name}.{region_name}.g.vcf.gz")
        wholegenome_combined_gvcf = os.path.join(cohort_dir, f"{cohort_name}.wholeGenome.g.vcf.gz")
        
        # Check if combined GVCF already exists in either format
        if os.path.exists(combined_gvcf):
            log(f"Combined GVCF already exists: {combined_gvcf}")
        elif os.path.exists(wholegenome_combined_gvcf):
            log(f"Combined GVCF already exists: {wholegenome_combined_gvcf}")
            combined_gvcf = wholegenome_combined_gvcf
        else:
            if len(cohort_gvcfs) == 1:
                sample, gvcf_path = cohort_gvcfs[0]
                shutil.copy(gvcf_path, combined_gvcf)
                if os.path.exists(gvcf_path + ".tbi"):
                    shutil.copy(gvcf_path + ".tbi", combined_gvcf + ".tbi")
                else:
                    # Index the combined GVCF if copy of index fails
                    index_vcf(combined_gvcf, args.threads)
                log(f"Only one sample in cohort {cohort_name}, copied GVCF directly to: {combined_gvcf}")
            else:
                # Combine GVCFs as usual
                combine_cmd = f"gatk CombineGVCFs -R {reference_fasta}"
                for sample, gvcf in cohort_gvcfs:
                    combine_cmd += f" --variant {gvcf}"
                combine_cmd += f" -O {combined_gvcf}"
                
                if not run(combine_cmd)[0]:
                    log(f"Failed to combine GVCFs for cohort {cohort_name}", "ERROR")
                    continue
                
                # Index the combined GVCF
                if not index_vcf(combined_gvcf, args.threads):
                    log(f"Failed to index combined GVCF for cohort {cohort_name}", "WARNING")
        
        # 2. Genotype GVCFs
        allsites_vcf = os.path.join(cohort_dir, f"{cohort_name}.allsites.geno.vcf.gz")
        if os.path.exists(allsites_vcf):
            log(f"Genotyped VCF already exists: {allsites_vcf}")
        else:
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
            
            # Index the genotyped VCF
            if not index_vcf(allsites_vcf, args.threads):
                log(f"Failed to index genotyped VCF for cohort {cohort_name}", "WARNING")
        
        # 3. Create SNP-only VCF
        final_filtered_vcf = os.path.join(cohort_dir, f"{cohort_name}.final.filtered.vcf.gz")
        if os.path.exists(final_filtered_vcf):
            log(f"Filtered VCF already exists: {final_filtered_vcf}")
        else:
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
            if not index_vcf(final_filtered_vcf, args.threads):
                log(f"Failed to index filtered VCF for cohort {cohort_name}", "ERROR")
                continue
        
        # 4. Create callable sites VCF
        allsites_filtered_vcf = os.path.join(cohort_dir, f"{cohort_name}.final.allsites.geno.filtered.vcf.gz")
        if os.path.exists(allsites_filtered_vcf):
            log(f"Callable sites VCF already exists: {allsites_filtered_vcf}")
        else:
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
            
            # Index the callable sites VCF
            if not index_vcf(allsites_filtered_vcf, args.threads):
                log(f"Failed to index callable sites VCF for cohort {cohort_name}", "WARNING")
        
        # 6. Split VCF into individual samples (moved before mask file creation to fix the error)
        per_sample_vcfs = []
        if os.path.exists(final_filtered_vcf):
            # Get sample names from VCF
            sample_list_cmd = f"bcftools query -l {final_filtered_vcf}"
            success, output = run(sample_list_cmd)
            if success and output:
                sample_names_in_vcf = output.strip().split()
                log(f"Sample names in VCF: {', '.join(sample_names_in_vcf)}")
                
                # Check if individual sample VCFs already exist
                for s in sample_names_in_vcf:
                    out_vcf = os.path.join(cohort_dir, f"{cohort_name}.{s}.vcf.gz")
                    if os.path.exists(out_vcf) and (os.path.exists(out_vcf + ".tbi") or os.path.exists(out_vcf + ".csi")):
                        log(f"Individual sample VCF already exists for {s}: {out_vcf}")
                        per_sample_vcfs.append(out_vcf)
                    else:
                        extract_cmd = f"bcftools view --samples {s} {final_filtered_vcf} --output-type z --output-file {out_vcf}"
                        if run(extract_cmd)[0]:
                            # Index the sample VCF
                            if index_vcf(out_vcf, args.threads):
                                per_sample_vcfs.append(out_vcf)
                                log(f"Created and indexed individual sample VCF for {s}")
                            else:
                                log(f"Failed to index individual sample VCF for {s}", "ERROR")
                        else:
                            log(f"Failed to create individual sample VCF for {s}", "ERROR")
            else:
                log(f"Failed to extract sample names from {final_filtered_vcf}", "ERROR")
        
        if not per_sample_vcfs:
            log(f"No sample VCFs available for cohort {cohort_name}", "ERROR")
            continue
        
        # 5. Create mask file
        mask_bed = os.path.join(cohort_dir, f"{cohort_name}.final.mask.bed")
        merged_mask_bed = os.path.join(cohort_dir, f"{cohort_name}.final.mask.merged.bed")
        final_mask = merged_mask_bed + ".gz"
        
        # Extract contigs from VCF before filtering mask
        vcf_chroms_tmp = os.path.join(cohort_dir, "vcf_chroms.txt")
        mask_filtered_tmp = merged_mask_bed + ".filtered"
        
        if os.path.exists(final_mask):
            log(f"Mask file already exists: {final_mask}")
        else:
            # Create mask BED
            if not os.path.exists(mask_bed):
                vcf2bed_cmd = f"zcat {allsites_filtered_vcf} | vcf2bed > {mask_bed}"
                if not run(vcf2bed_cmd)[0]:
                    log(f"Failed to convert VCF to BED for cohort {cohort_name}", "ERROR")
                    continue
            
            sorted_mask_bed = mask_bed.replace(".bed", ".sorted.bed")
            sort_cmd = f"sort -k1,1 -k2,2n {mask_bed} > {sorted_mask_bed}"
            if not run(sort_cmd)[0]:
                log(f"Failed to sort BED file for cohort {cohort_name}", "ERROR")
                continue

            # Merge BED intervals
            if not os.path.exists(merged_mask_bed):
                merge_cmd = f"bedtools merge -i {sorted_mask_bed} > {merged_mask_bed}"
                if not run(merge_cmd)[0]:
                    log(f"Failed to merge BED intervals for cohort {cohort_name}", "ERROR")
                    continue

                # Ensure merged BED is properly sorted before compression (fix for AssertionError)
                resort_merged = merged_mask_bed.replace(".bed", ".resorted.bed")
                sort_resort_cmd = f"sort -k1,1 -k2,2n {merged_mask_bed} > {resort_merged}"
                if not run(sort_resort_cmd)[0]:
                    log(f"Failed to resort merged BED for cohort {cohort_name}", "ERROR")
                    continue
                shutil.move(resort_merged, merged_mask_bed)
            
            # Filter mask to include only contigs in VCF
            log(f"Filtering mask to match VCF contigs for cohort {cohort_name}")
            
            # Extract contigs from VCF
            vcf_chroms_cmd = f"bcftools query -f '%CHROM\\n' {per_sample_vcfs[0]} | sort | uniq > {vcf_chroms_tmp}"
            if not run(vcf_chroms_cmd)[0]:
                log("Failed to extract contigs from VCF", "ERROR")
                continue

            # Filter the merged BED file to only include those contigs
            filter_mask_cmd = (
                f"grep -Ff {vcf_chroms_tmp} {merged_mask_bed} > {mask_filtered_tmp}"
            )
            if not run(filter_mask_cmd)[0]:
                log("Failed to filter mask BED to VCF contigs", "ERROR")
                continue

            # Replace original merged mask with filtered version
            shutil.move(mask_filtered_tmp, merged_mask_bed)
            log("Filtered mask BED now matches VCF contigs")
            
            # Final resort to ensure proper ordering within each contig
            final_resort = merged_mask_bed + ".final_sort"
            final_sort_cmd = f"sort -k1,1 -k2,2n {merged_mask_bed} > {final_resort}"
            if not run(final_sort_cmd)[0]:
                log("Failed to do final resort of BED file", "ERROR")
                continue
            shutil.move(final_resort, merged_mask_bed)
            
            # Compress mask
            if not run(f"gzip -f {merged_mask_bed}")[0]:
                log(f"Failed to compress mask BED for cohort {cohort_name}", "ERROR")
                continue
        
        # 7. Generate Multihetsep file
        multihetsep_suffix = region_name
        multihetsep_out = os.path.join(cohort_dir, f"{cohort_name}.{multihetsep_suffix}.mhs")
        wholegenome_multihetsep_out = os.path.join(cohort_dir, f"{cohort_name}.wholeGenome.mhs")
        
        # Check if Multihetsep file already exists in either format
        if os.path.exists(multihetsep_out):
            log(f"Multihetsep file already exists: {multihetsep_out}")
        elif os.path.exists(wholegenome_multihetsep_out):
            log(f"Multihetsep file already exists: {wholegenome_multihetsep_out}")
            multihetsep_out = wholegenome_multihetsep_out
        else:
            # Download generate_multihetsep.py script if needed
            gen_mhs_script = os.path.join(cohort_dir, "generate_multihetsep.py")
            if not os.path.exists(gen_mhs_script):
                log("Downloading generate_multihetsep.py script")
                wget_cmd = f"wget -O {gen_mhs_script} https://raw.githubusercontent.com/stschiff/msmc-tools/master/generate_multihetsep.py"
                if not run(wget_cmd)[0] or not run(f"chmod u+x {gen_mhs_script}")[0]:
                    log("Failed to download or set permissions for generate_multihetsep.py", "ERROR")
                    continue
            
            # Patch the script to handle unsorted positions
            if not patch_generate_multihetsep(gen_mhs_script):
                log("Failed to patch generate_multihetsep.py script", "ERROR")
                continue
            
            # If skipping chromosome splitting, process entire VCF at once
            if args.skip_chromosome_splitting:
                log("Processing all chromosomes at once (skipping chromosome splitting)")
                
                # Generate multihetsep for the whole genome
                mhs_cmd = f"python3 {gen_mhs_script} --mask={final_mask} " + " ".join(per_sample_vcfs) + f" > {multihetsep_out}"
                if run(mhs_cmd)[0]:
                    log(f"Successfully generated multihetsep file: {multihetsep_out}")
                else:
                    log("Failed to generate multihetsep file", "ERROR")
                    continue
            else:
                # Process one chromosome at a time
                # Get all chromosomes/contigs from VCF
                chroms_cmd = f"bcftools query -f '%CHROM\\n' {per_sample_vcfs[0]} | sort | uniq"
                success, chroms_output = run(chroms_cmd)
                if not success or not chroms_output:
                    log("Failed to get chromosomes from VCF", "ERROR")
                    continue
                
                chromosomes = chroms_output.strip().split('\n')
                log(f"Processing {len(chromosomes)} chromosomes/contigs")
                
                # Create a temporary directory for chromosome-specific files
                tmp_dir = os.path.join(cohort_dir, "tmp")
                os.makedirs(tmp_dir, exist_ok=True)
                
                # Process each chromosome separately
                chrom_mhs_files = []
                for i, chrom in enumerate(chromosomes):
                    log(f"Processing chromosome/contig {i+1}/{len(chromosomes)}: {chrom}")
                    chrom_mhs = os.path.join(tmp_dir, f"{chrom}.mhs")
                    
                    # Extract chromosome region from VCF
                    chrom_vcfs = []
                    for vcf in per_sample_vcfs:
                        sample = os.path.basename(vcf).split('.')[-2]
                        chrom_vcf = os.path.join(tmp_dir, f"{sample}.{chrom}.vcf.gz")
                        
                        # Check if VCF is indexed before trying to extract
                        if not (os.path.exists(vcf + ".tbi") or os.path.exists(vcf + ".csi")):
                            log(f"VCF is not indexed, creating index for {vcf}", "WARNING")
                            if not index_vcf(vcf, args.threads):
                                log(f"Failed to index VCF {vcf}, skipping chromosome {chrom}", "ERROR")
                                continue
                                
                        # Now extract chromosome
                        extract_cmd = f"bcftools view -r {chrom} {vcf} -Oz -o {chrom_vcf}"
                        if run(extract_cmd)[0]:
                            # Make sure chromosome VCF is indexed
                            if index_vcf(chrom_vcf, 1):
                                chrom_vcfs.append(chrom_vcf)
                            else:
                                log(f"Failed to index chromosome VCF for {chrom}", "ERROR")
                        else:
                            log(f"Failed to extract chromosome {chrom} from VCF", "WARNING")
                    
                    if not chrom_vcfs:
                        log(f"No VCFs available for chromosome {chrom}, skipping", "WARNING")
                        continue
                    
                    # Extract chromosome region from mask
                    chrom_mask = os.path.join(tmp_dir, f"{chrom}.mask.bed.gz")
                    extract_mask_cmd = f"gunzip -c {final_mask} | grep -P '^{chrom}\\t' | gzip > {chrom_mask}"
                    if not run(extract_mask_cmd)[0]:
                        log(f"Failed to extract chromosome {chrom} from mask", "WARNING")
                        # Try without a mask
                        mhs_cmd = f"python3 {gen_mhs_script} " + " ".join(chrom_vcfs) + f" > {chrom_mhs}"
                    else:
                        # Generate multihetsep for this chromosome
                        mhs_cmd = f"python3 {gen_mhs_script} --mask={chrom_mask} " + " ".join(chrom_vcfs) + f" > {chrom_mhs}"
                    
                    if run(mhs_cmd)[0]:
                        # Check if file is empty
                        if os.path.getsize(chrom_mhs) > 0:
                            chrom_mhs_files.append(chrom_mhs)
                        else:
                            log(f"Generated empty multihetsep file for chromosome {chrom}", "WARNING")
                    else:
                        log(f"Failed to generate multihetsep for chromosome {chrom}", "WARNING")
                
                # Combine all chromosome multihetsep files
                if chrom_mhs_files:
                    concat_cmd = f"cat {' '.join(chrom_mhs_files)} > {multihetsep_out}"
                    if not run(concat_cmd)[0]:
                        log("Failed to concatenate chromosome multihetsep files", "ERROR")
                    else:
                        log(f"Successfully generated multihetsep file: {multihetsep_out}")
                    
                    # Clean up temporary files
                    if run(f"rm -rf {tmp_dir}")[0]:
                        log("Cleaned up temporary files")
                else:
                    # If processing by chromosome failed, try processing whole VCF at once as fallback
                    log("No chromosome multihetsep files were generated, trying without chromosome splitting", "WARNING")
                    mhs_cmd = f"python3 {gen_mhs_script} --mask={final_mask} " + " ".join(per_sample_vcfs) + f" > {multihetsep_out}"
                    if run(mhs_cmd)[0]:
                        log(f"Successfully generated multihetsep file using whole-VCF method: {multihetsep_out}")
                    else:
                        log("Failed to generate multihetsep file", "ERROR")
                        continue
        
        # 8. Create subset if requested
        if args.create_subset and os.path.exists(multihetsep_out):
            subset_out = os.path.join(
                cohort_dir, 
                f"{cohort_name}.{multihetsep_suffix}.{args.subset_start/1000000:.0f}Mbto{args.subset_end/1000000:.0f}Mb.subset.mhs"
            )
            if os.path.exists(subset_out):
                log(f"Subset file already exists: {subset_out}")
            else:
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
        if os.path.exists(multihetsep_out):
            check_cmd = f"head -n 10 {multihetsep_out}"
            run(check_cmd)
        else:
            log(f"Cannot perform quality check - Multihetsep file not found", "WARNING")  # noqa: F541
        
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