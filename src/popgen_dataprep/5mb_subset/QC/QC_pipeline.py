import os
import subprocess
import logging
from Bio import SeqIO
import matplotlib.pyplot as plt
import pandas as pd

# Set up logging
def setup_logging(log_directory, log_filename="pipeline.log"):
    """Configure logging to write to both file and console with detailed formatting."""
    os.makedirs(log_directory, exist_ok=True)
    log_file_path = os.path.join(log_directory, log_filename)
    
    # Create logger
    logger = logging.getLogger('bioinfo_pipeline')
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

# Number of CPU cores for mapping
NUM_THREADS = 30

# Define paths
BASE_DIR = "/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/Popgen_analysis/dataprep/5mb_subset"
QC_DIR = os.path.join(BASE_DIR, "QC")
INTERMEDIATE_DIR = os.path.join(QC_DIR, "Intermediate_data")

# Ensure top-level directories exist
os.makedirs(QC_DIR, exist_ok=True)
os.makedirs(INTERMEDIATE_DIR, exist_ok=True)

# Initialize logging
logger, log_file_path = setup_logging(QC_DIR)
logger.info(f"🔍 Logging initialized. Log file: {log_file_path}")
logger.info(f"🛠️ Pipeline started with {NUM_THREADS} threads")
logger.info(f"📂 Base directory: {BASE_DIR}")

# Define input files
samples = {
    "HT": {  # Long-read Nanopore (ONT dataset candidate)
        "assembly": os.path.join(BASE_DIR, "HT_assembly_5mb_subset.fasta"),
        "reads": "/data/proj2/home/students/m.borgmann/Master_thesis/data/raw/Reference_data/WGS_projects/Hexaplex_trunculus_nanopore/20241217_DNA_Tellier_HT_31_PB/20241217_1155_2A_PAY72199_a0bfda36/fastq_pass/*.fastq.gz",
        "mapper": "minimap2",
        "data_type": "ONT"
    },
    "HT2": {  # Short-read Illumina
        "assembly": os.path.join(BASE_DIR, "HT2_assembly_5mb_subset.fasta"),
        "reads_1": "/data/proj2/home/students/m.borgmann/Master_thesis/data/raw/Reference_data/WGS_projects/PRJNA1106542_Hexaplex_trunculus/SRR28865916/SRR28865916_R1.fastq",
        "reads_2": "/data/proj2/home/students/m.borgmann/Master_thesis/data/raw/Reference_data/WGS_projects/PRJNA1106542_Hexaplex_trunculus/SRR28865916/SRR28865916_R2.fastq",
        "mapper": "bwa-mem2",
        "data_type": "ILLUMINA"
    },
    "BB": {  # Short-read Illumina
        "assembly": os.path.join(BASE_DIR, "BB_assembly_5mb_subset.fasta"),
        "reads_1": "/data/proj2/home/students/m.borgmann/Master_thesis/data/raw/Reference_data/WGS_projects/PRJNA1106534_Bolinus_brandaris/SRR28863561/SRR28863561_R1.fastq",
        "reads_2": "/data/proj2/home/students/m.borgmann/Master_thesis/data/raw/Reference_data/WGS_projects/PRJNA1106534_Bolinus_brandaris/SRR28863561/SRR28863561_R2.fastq",
        "mapper": "bwa-mem2",
        "data_type": "ILLUMINA"
    },
}

# List available datasets and prompt for ONT datasets
sample_keys = list(samples.keys())
logger.info("Available datasets:")
for idx, key in enumerate(sample_keys, start=1):
    logger.info(f"{idx}: {key}")

ont_input = input("Enter dataset numbers (semicolon separated) that are ONT, or press Enter if none: ")

# Build a set of dataset keys that are ONT based on user input.
ont_datasets = set()
if ont_input.strip():
    numbers = [num.strip() for num in ont_input.split(";") if num.strip().isdigit()]
    for num in numbers:
        idx = int(num) - 1  # convert to 0-based index
        if 0 <= idx < len(sample_keys):
            ont_datasets.add(sample_keys[idx])
logger.info(f"Datasets marked as ONT: {ont_datasets}")

####################################################
# Create per-sample intermediate directories
####################################################
sample_intermediate_dirs = {}
for sample in samples:
    sample_dir = os.path.join(INTERMEDIATE_DIR, sample)
    os.makedirs(sample_dir, exist_ok=True)
    sample_intermediate_dirs[sample] = sample_dir
    logger.debug(f"Created intermediate directory for {sample}: {sample_dir}")

# Paths for raw VCF (unfiltered) and final filtered VCF, BAM and depth files.
raw_bam_files = {}        
filtered_bam_files = {}   
bam_files = {}
depth_files = {}
raw_vcf_files = {}
filtered_vcf_files = {}

for sample in samples:
    sample_dir = sample_intermediate_dirs[sample]
    
    # We'll first produce a "raw" BAM from alignment, then a "filtered" BAM removing secondaries, etc.
    raw_bam_files[sample] = os.path.join(sample_dir, f"{sample}_raw.bam")        
    filtered_bam_files[sample] = os.path.join(sample_dir, f"{sample}_filtered.bam")  
    bam_files[sample] = os.path.join(sample_dir, f"{sample}_sorted.bam")
    depth_files[sample] = os.path.join(sample_dir, f"{sample}_depth.txt")
    raw_vcf_files[sample] = os.path.join(sample_dir, f"{sample}_raw.vcf.gz")
    filtered_vcf_files[sample] = os.path.join(sample_dir, f"{sample}_filtered.vcf.gz")
    logger.debug(f"File paths for {sample} defined")

qc_report = os.path.join(QC_DIR, "QC_report.txt")
logger.info(f"QC report will be saved to: {qc_report}")

####################################################
# Utility Functions
####################################################
def median(lst):
    """Calculate median of a list."""
    n = len(lst)
    if n == 0:
        return 0
    s = sorted(lst)
    mid = n // 2
    if n % 2 == 1:
        return s[mid]
    else:
        return (s[mid - 1] + s[mid]) / 2

def compute_heterozygosity_percent(vcf_path, het_count):
    """
    Fraction of sites that are heterozygous = (# HET / total # called sites) * 100.
    Because we use 'bcftools call -m' (not -v), the VCF should include ref calls too.
    """
    try:
        total_sites = int(subprocess.check_output(
            f"bcftools view -H {vcf_path} | wc -l", shell=True
        ).strip())
        logger.debug(f"Total sites in {vcf_path}: {total_sites}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to count total sites in {vcf_path}: {e}")
        total_sites = 0
    if total_sites > 0:
        het_percent = (het_count / total_sites) * 100
        logger.debug(f"Heterozygosity calculation: {het_count}/{total_sites} = {het_percent:.2f}%")
        return het_percent
    else:
        logger.warning(f"No sites found in {vcf_path}, returning 0% heterozygosity")
        return 0

def cleanup_bwa_mem2_files(assembly_path):
    """Remove intermediary index files created by bwa-mem2 indexing."""
    logger.info(f"Cleaning up bwa-mem2 index files for {assembly_path}")
    suffixes = [".0123", ".amb", ".ann", ".bwt.2bit.64", ".fai", ".pac"]
    for suf in suffixes:
        file_to_remove = assembly_path + suf
        if os.path.exists(file_to_remove):
            os.remove(file_to_remove)
            logger.debug(f"Removed {file_to_remove}")
        else:
            logger.debug(f"Index file not found: {file_to_remove}")

def run(cmd, workdir=None):
    """Execute a shell command and log the output."""
    logger.info(f"🔧 Running: {cmd}")
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
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed with exit code {e.returncode}")
        logger.error(f"Command stdout: {e.stdout.strip() if e.stdout else 'None'}")
        logger.error(f"Command stderr: {e.stderr.strip() if e.stderr else 'None'}")
        return False

####################################################
# Step 1A: Map Reads -> raw.bam
####################################################
def map_reads(sample, assembly, raw_bam, mapper):
    """Map reads using minimap2 (ONT) or bwa-mem2 (Illumina). Output = raw BAM."""
    logger.info(f"🚀 Mapping reads for {sample} using {mapper} with {NUM_THREADS} threads...")
    sam_file = raw_bam.replace(".bam", ".sam")
    logger.debug(f"Temporary SAM file: {sam_file}")

    if mapper == "minimap2":
        # ONT
        logger.debug(f"Using minimap2 for ONT reads: {samples[sample]['reads']}")
        cmd = f"minimap2 -t {NUM_THREADS} -ax map-ont {assembly} {samples[sample]['reads']} > {sam_file}"
    elif mapper == "bwa-mem2":
        # Illumina
        logger.debug(f"Using bwa-mem2 for Illumina reads: R1={samples[sample]['reads_1']}, R2={samples[sample]['reads_2']}")
        cmd = (
            f"bwa-mem2 index {assembly} && "
            f"bwa-mem2 mem -t {NUM_THREADS} {assembly} "
            f"{samples[sample]['reads_1']} {samples[sample]['reads_2']} > {sam_file}"
        )
    else:
        error_msg = f"Unknown mapper specified for {sample}: {mapper}"
        logger.error(error_msg)
        raise ValueError(error_msg)

    # Run the alignment -> SAM
    success = run(cmd)
    if not success:
        logger.error(f"Alignment failed for {sample}")
        return None

    # Convert SAM to an unsorted BAM
    logger.info(f"Converting SAM to BAM for {sample}")
    bam_cmd = f"samtools view -@ {NUM_THREADS} -bS {sam_file} -o {raw_bam}"
    success = run(bam_cmd)
    if not success:
        logger.error(f"SAM to BAM conversion failed for {sample}")
        return None

    # Remove intermediate SAM file
    logger.debug(f"Removing temporary SAM file: {sam_file}")
    os.remove(sam_file)
    logger.info(f"✅ Read mapping completed for {sample}")

    return raw_bam

####################################################
# Step 1B: Filter raw.bam -> filtered.bam
####################################################
def filter_bam(sample, raw_bam, filtered_bam):
    """
    Filter out undesired alignments:
      - ONT: remove secondary(0x100) + supplementary(0x800)
      - Illumina: optionally remove duplicates via samtools markdup
    """
    data_type = samples[sample]["data_type"]
    temp_bam = filtered_bam.replace(".bam", ".tmp.bam")
    logger.info(f"🔍 Filtering BAM file for {sample} ({data_type} data type)")

    if data_type == "ONT":
        # Filter out secondary + supplementary
        # -F 0x900 means: skip reads with either 0x100 or 0x800 bits set
        logger.debug(f"Applying ONT-specific filtering: removing secondary/supplementary alignments")
        cmd_filter = f"samtools view -@ {NUM_THREADS} -F 0x900 -b {raw_bam} -o {filtered_bam}"
        if not run(cmd_filter):
            logger.error(f"ONT filtering failed for {sample}")
            return None

    elif data_type == "ILLUMINA":
        # (1) Optionally mark duplicates
        logger.debug(f"Applying Illumina-specific processing: marking duplicates")
        
        # 1a. Sort by name for fixmates
        logger.debug(f"Name-sorting BAM for fixmate")
        cmd_sort_name = f"samtools sort -@ {NUM_THREADS} -n {raw_bam} -o {temp_bam}"
        if not run(cmd_sort_name):
            logger.error(f"Name sorting failed for {sample}")
            return None

        logger.debug(f"Running fixmate")
        cmd_fixmate = f"samtools fixmate -@ {NUM_THREADS} -m {temp_bam} {temp_bam}.fixmate.bam"
        if not run(cmd_fixmate):
            logger.error(f"Fixmate failed for {sample}")
            return None
        
        # 1c. Sort by coordinate
        logger.debug(f"Coordinate-sorting BAM")
        cmd_sort_coord = f"samtools sort -@ {NUM_THREADS} {temp_bam}.fixmate.bam -o {temp_bam}.sorted.bam"
        if not run(cmd_sort_coord):
            logger.error(f"Coordinate sorting failed for {sample}")
            return None
        
        # 1d. Mark duplicates
        logger.debug(f"Marking duplicates")
        cmd_markdup = f"samtools markdup -@ {NUM_THREADS} {temp_bam}.sorted.bam {temp_bam}.markdup.bam"
        if not run(cmd_markdup):
            logger.error(f"Marking duplicates failed for {sample}")
            return None

        # 2. Remove duplicates ( -F 0x400 ) and secondary/supp ( -F 0x900 )
        logger.debug(f"Filtering out duplicates and secondary/supplementary alignments")
        cmd_view = (
            "samtools view -@ {threads} -b -F 0x400 -F 0x900 {inbam} -o {outbam}"
            .format(
                threads=NUM_THREADS,
                inbam=f"{temp_bam}.markdup.bam",
                outbam=filtered_bam
            )
        )
        if not run(cmd_view):
            logger.error(f"Filtering duplicates failed for {sample}")
            return None

        # Cleanup
        logger.debug(f"Cleaning up temporary files")
        for f in [temp_bam, f"{temp_bam}.fixmate.bam", f"{temp_bam}.sorted.bam", f"{temp_bam}.markdup.bam"]:
            if os.path.exists(f):
                os.remove(f)
                logger.debug(f"Removed {f}")

    else:
        # Unrecognized data type, just copy raw -> filtered
        logger.warning(f"Unrecognized data type '{data_type}' for {sample}, just copying raw BAM to filtered")
        cmd_cp = f"cp {raw_bam} {filtered_bam}"
        if not run(cmd_cp):
            logger.error(f"Copying BAM failed for {sample}")
            return None

    logger.info(f"✅ BAM filtering completed for {sample}")
    return filtered_bam

####################################################
# Step 1C: Sort + Index filtered.bam -> final .bam
####################################################
def sort_and_index_bam(sample, filtered_bam, final_bam):
    """Sort and index the filtered BAM for coverage & variant calling."""
    logger.info(f"📚 Sorting & indexing filtered BAM for {sample}")
    
    # Sort BAM
    logger.debug(f"Sorting BAM file")
    cmd_sort = f"samtools sort -@ {NUM_THREADS} {filtered_bam} -o {final_bam}"
    if not run(cmd_sort):
        logger.error(f"BAM sorting failed for {sample}")
        return None
    
    # Index BAM
    logger.debug(f"Indexing BAM file")
    cmd_index = f"samtools index -@ {NUM_THREADS} {final_bam}"
    if not run(cmd_index):
        logger.error(f"BAM indexing failed for {sample}")
        return None
    
    logger.info(f"✅ BAM sorting and indexing completed for {sample}")
    return final_bam

####################################################
# Step 2: Compute Coverage (raw), then truncate in Python
####################################################
def compute_coverage(sample, bam, depth_output):
    """
    Compute coverage statistics without any capping by samtools.
    We'll do the capping/truncation in Python for the QC report.
    """
    logger.info(f"📊 Computing coverage statistics for {sample} (raw coverage)...")

    # Run samtools depth WITHOUT -d so we capture outliers
    cmd_depth = f"samtools depth -@ {NUM_THREADS} {bam} > {depth_output}"
    if not run(cmd_depth):
        logger.error(f"Depth calculation failed for {sample}")
        return 0, 0, 0

    # Load raw coverage from output
    logger.debug(f"Processing raw depth file: {depth_output}")
    try:
        with open(depth_output, "r") as f:
            depths = [int(line.split()[2]) for line in f if len(line.split()) >= 3]
    except Exception as e:
        logger.error(f"Error reading depth file for {sample}: {e}")
        return 0, 0, 0

    if not depths:
        logger.warning(f"No depth data found for {sample}")
        return 0, 0, 0

    # Decide coverage cap for truncation in QC stats (Illumina=250, ONT=500)
    cap = 250 if samples[sample]['data_type'] == 'ILLUMINA' else 500
    truncated_depths = [min(d, cap) for d in depths]

    # Compute summary stats using truncated depths
    mean_cov = sum(truncated_depths) / len(truncated_depths)
    med_cov = median(truncated_depths)
    ten_x_sites = (sum(1 for d in truncated_depths if d >= 10) / len(truncated_depths)) * 100

    logger.info(f"Coverage statistics for {sample} (using truncated depths at {cap}x for QC):")
    logger.info(f"  - Mean coverage: {mean_cov:.2f}x")
    logger.info(f"  - Median coverage: {med_cov:.2f}x")
    logger.info(f"  - Sites with ≥10x coverage: {ten_x_sites:.2f}%")
    logger.debug(f"  - Total positions analyzed: {len(depths)}")

    return mean_cov, med_cov, ten_x_sites

####################################################
# Step 3: Call Variants with BCFtools (raw VCF)
####################################################
def call_variants(sample, assembly, bam, vcf_output):
    """
    Call variants using bcftools with per-genotype fields (AD, DP, SP)
    and include non-variant (reference) sites in the VCF.
    
    For ONT data, coverage is capped at 500; for Illumina, it's capped at 250
    *here in bcftools mpileup* so as not to blow up read counts in variant calling.
    """
    logger.info(f"🧬 Calling variants for {sample}...")

    # Define mpileup options based on data type.
    if sample in ont_datasets:
        mpileup_options = (
            f"-X ont -d 500 -Ou -f {assembly} {bam} "  # ONT-specific: cap depth at 500
            f"-a AD,DP,SP "                           # Genotype annotations
            f"--threads {NUM_THREADS}"
        )
    else:
        mpileup_options = (
            f"-d 250 -Ou -f {assembly} {bam} "         # Illumina-specific: cap depth at 250
            f"-a AD,DP,SP "                           # Genotype annotations
            f"--threads {NUM_THREADS}"
        )

    # bcftools pipeline
    mpileup_cmd = f"bcftools mpileup {mpileup_options}"
    call_cmd = f"bcftools call --ploidy 2 -m -A -Oz -o {vcf_output}"
    full_cmd = f"{mpileup_cmd} | {call_cmd}"

    logger.debug(f"Running variant calling pipeline: {full_cmd}")
    success = run(full_cmd)
    if not success:
        logger.error(f"Variant calling pipeline failed for {sample}")
        return False

    logger.debug("Indexing VCF file")
    index_success = run(f"bcftools index --threads {NUM_THREADS} {vcf_output}")
    if not index_success:
        logger.error(f"VCF indexing failed for {sample}")
        return False

    logger.info(f"✅ Variant calling completed for {sample}")
    return True

####################################################
# Step 4: Filter Variants (DP>=10, QUAL>=30)
####################################################
def filter_variants(sample, raw_vcf, filtered_vcf):
    """
    Filter variants with DP >= 10 and QUAL >= 30.
    Because this is a single-sample VCF, bcftools places DP in INFO for each site by default.
    """
    logger.info(f"🚦 Filtering variants for {sample} with DP>=10 and QUAL>=30...")
    
    filter_expr = "'QUAL<30 || INFO/DP<10'"
    filter_cmd = f"bcftools filter -e {filter_expr} {raw_vcf} -Oz -o {filtered_vcf}"
    
    try:
        logger.debug(f"Running bcftools filter with expression: {filter_expr}")
        success = run(filter_cmd)
        if not success:
            logger.error(f"Variant filtering failed for {sample}")
            return False
            
        logger.debug(f"Indexing filtered VCF")
        index_success = run(f"bcftools index --threads {NUM_THREADS} {filtered_vcf}")
        if not index_success:
            logger.error(f"Filtered VCF indexing failed for {sample}")
            return False
            
        logger.info(f"✅ Variant filtering completed for {sample}")
        return True
    except Exception as e:
        logger.error(f"❌ Unexpected error in variant filtering for {sample}: {e}")
        return False

####################################################
# Helper function: Count Heterozygous & Missing Sites
####################################################
def count_het_and_missing(vcf_path):
    """
    Count heterozygous sites and missing sites in a (filtered) VCF.
    Because we used 'bcftools call -m', these are among ALL called sites (ref + alt).
    """
    if not os.path.exists(vcf_path):
        logger.warning(f"VCF file not found: {vcf_path}")
        return 0, 0

    logger.debug(f"Counting heterozygous and missing sites in {vcf_path}")
    
    # Number of lines with a heterozygous genotype
    try:
        het_count = int(subprocess.check_output(
            f"bcftools view -g het {vcf_path} | wc -l", shell=True
        ).strip())
        logger.debug(f"Heterozygous sites: {het_count}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to count heterozygous sites: {e}")
        het_count = 0

    # Number of lines with a missing genotype './.'
    try:
        missing_count = int(subprocess.check_output(
            f"bcftools query -f '%CHROM\\t%POS\\t[%GT]\\n' {vcf_path} | grep '\\./\\.' | wc -l",
            shell=True
        ).strip())
        logger.debug(f"Missing sites: {missing_count}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to count missing sites: {e}")
        missing_count = 0

    return het_count, missing_count

####################################################
# Run the Pipeline for Each Sample
####################################################
logger.info("\n=============================================")
logger.info("🚀 Starting pipeline execution for all samples")
logger.info("=============================================\n")

results = {}
for sample in samples:
    logger.info(f"\n=======================")
    logger.info(f"🚀 Processing sample: {sample}")
    logger.info(f"=======================\n")

    assembly = samples[sample]["assembly"]
    mapper = samples[sample]["mapper"]
    logger.debug(f"Sample details - Assembly: {assembly}, Mapper: {mapper}, Data type: {samples[sample]['data_type']}")

    # 1A) Map Reads -> raw.bam
    if os.path.exists(raw_bam_files[sample]):
        logger.info(f"✅ Raw BAM found for {sample}. Skipping mapping.")
        raw_bam = raw_bam_files[sample]
    else:
        logger.info(f"Raw BAM not found. Starting read mapping for {sample}.")
        raw_bam = map_reads(sample, assembly, raw_bam_files[sample], mapper)
        if not raw_bam:
            logger.error(f"❌ Read mapping failed for {sample}. Skipping further processing.")
            results[sample] = {
                "Mean Coverage": 0,
                "Median Coverage": 0,
                "% Sites ≥10x Coverage": 0,
                "Heterozygous Sites": 0,
                "Heterozygosity (%)": 0,
                "Missing Sites": 0,
                "Status": "Failed - Mapping error"
            }
            continue

    # 1B) Filter raw.bam -> filtered.bam
    if os.path.exists(filtered_bam_files[sample]):
        logger.info(f"✅ Filtered BAM found for {sample}. Skipping filtering step.")
        filtered_bam = filtered_bam_files[sample]
    else:
        logger.info(f"Filtered BAM not found. Starting BAM filtering for {sample}.")
        filtered_bam = filter_bam(sample, raw_bam, filtered_bam_files[sample])
        if not filtered_bam:
            logger.error(f"❌ BAM filtering failed for {sample}. Skipping further processing.")
            results[sample] = {
                "Mean Coverage": 0,
                "Median Coverage": 0,
                "% Sites ≥10x Coverage": 0,
                "Heterozygous Sites": 0,
                "Heterozygosity (%)": 0,
                "Missing Sites": 0,
                "Status": "Failed - Filtering error"
            }
            continue

    # 1C) Sort + Index -> final .bam
    if os.path.exists(bam_files[sample]) and os.path.exists(bam_files[sample] + ".bai"):
        logger.info(f"✅ Final sorted BAM found for {sample}. Skipping sort/index.")
        final_bam = bam_files[sample]
    else:
        logger.info(f"Sorted BAM not found. Starting BAM sorting for {sample}.")
        final_bam = sort_and_index_bam(sample, filtered_bam, bam_files[sample])
        if not final_bam:
            logger.error(f"❌ BAM sorting/indexing failed for {sample}. Skipping further processing.")
            results[sample] = {
                "Mean Coverage": 0,
                "Median Coverage": 0,
                "% Sites ≥10x Coverage": 0,
                "Heterozygous Sites": 0,
                "Heterozygosity (%)": 0,
                "Missing Sites": 0,
                "Status": "Failed - Sorting error"
            }
            continue

    # 2) Compute Coverage
    if os.path.exists(depth_files[sample]):
        logger.info(f"✅ Depth file found for {sample}. Loading and re‐computing truncated coverage data.")
        try:
            with open(depth_files[sample], "r") as f:
                raw_depths = [int(line.split()[2]) for line in f if len(line.split()) >= 3]
        except Exception as e:
            logger.error(f"Error reading existing depth file for {sample}: {e}")
            logger.info("Recomputing coverage from BAM.")
            mean_cov, median_cov, ten_x_sites = compute_coverage(sample, final_bam, depth_files[sample])
        else:
            if raw_depths:
                # Truncate coverage for QC stats
                cap = 250 if samples[sample]['data_type'] == 'ILLUMINA' else 500
                truncated = [min(d, cap) for d in raw_depths]
                mean_cov = sum(truncated)/len(truncated)
                median_cov = median(truncated)
                ten_x_sites = (sum(1 for d in truncated if d >= 10)/len(truncated))*100
                logger.info(f"Loaded coverage statistics for {sample} (truncated at {cap}x):")
                logger.info(f"  - Mean coverage: {mean_cov:.2f}x")
                logger.info(f"  - Median coverage: {median_cov:.2f}x")
                logger.info(f"  - Sites with ≥10x coverage: {ten_x_sites:.2f}%")
            else:
                logger.warning(f"Depth file exists but contains no valid data for {sample}")
                mean_cov, median_cov, ten_x_sites = 0, 0, 0
    else:
        logger.info(f"Depth file not found. Computing coverage for {sample}.")
        mean_cov, median_cov, ten_x_sites = compute_coverage(sample, final_bam, depth_files[sample])

    # 3) Variant Calling
    raw_vcf = raw_vcf_files[sample]
    if os.path.exists(raw_vcf):
        logger.info(f"✅ Raw VCF found for {sample}. Skipping bcftools call.")
    else:
        logger.info(f"Raw VCF not found. Calling variants for {sample}.")
        success = call_variants(sample, assembly, final_bam, raw_vcf)
        if not success:
            logger.error(f"❌ Variant calling step failed for {sample}.")
            results[sample] = {
                "Mean Coverage": mean_cov,
                "Median Coverage": median_cov,
                "% Sites ≥10x Coverage": ten_x_sites,
                "Heterozygous Sites": 0,
                "Heterozygosity (%)": 0,
                "Missing Sites": 0,
                "Status": "Failed - Variant calling error"
            }
            continue

    # 4) Filter VCF
    filtered_vcf = filtered_vcf_files[sample]
    if os.path.exists(filtered_vcf):
        logger.info(f"✅ Filtered VCF found for {sample}. Skipping filtering.")
    else:
        logger.info(f"Filtered VCF not found. Filtering variants for {sample}.")
        success = filter_variants(sample, raw_vcf, filtered_vcf)
        if not success:
            logger.error(f"❌ Filtering step failed for {sample}.")
            results[sample] = {
                "Mean Coverage": mean_cov,
                "Median Coverage": median_cov,
                "% Sites ≥10x Coverage": ten_x_sites,
                "Heterozygous Sites": 0,
                "Heterozygosity (%)": 0,
                "Missing Sites": 0,
                "Status": "Failed - VCF filtering error"
            }
            continue

    # 5) Count heterozygous & missing sites
    het_count, missing_count = count_het_and_missing(filtered_vcf)
    het_percent = compute_heterozygosity_percent(filtered_vcf, het_count)

    results[sample] = {
        "Mean Coverage": mean_cov,
        "Median Coverage": median_cov,
        "% Sites ≥10x Coverage": ten_x_sites,
        "Heterozygous Sites": het_count,
        "Heterozygosity (%)": het_percent,
        "Missing Sites": missing_count,
        "Status": "Completed successfully"
    }
    logger.info(f"✅ Processing completed for sample {sample}")

    # 6) Cleanup BWA mem2 intermediary files for sample BB only (optional)
    if sample == "BB" and mapper == "bwa-mem2":
        logger.info("🧹 Cleaning up intermediary bwa-mem2 files for sample BB...")
        cleanup_bwa_mem2_files(assembly)

####################################################
# Step 5: Save QC Report with a neat table
####################################################
logger.info("\n=======================")
logger.info("📊 Generating QC report")
logger.info("=======================\n")

header = (
    "Sample".ljust(10) +
    "Mean Cov".rjust(10) +
    "Median Cov".rjust(12) +
    "% ≥10x".rjust(10) +
    "Het Sites".rjust(12) +
    "Het (%)".rjust(10) +
    "Missing".rjust(10) +
    "Status".rjust(20)
)
separator = "-" * len(header)
report_lines = [header, separator]
for sample, stats in results.items():
    line = (
        sample.ljust(10) +
        f"{stats['Mean Coverage']:.2f}".rjust(10) +
        f"{stats['Median Coverage']:.2f}".rjust(12) +
        f"{stats['% Sites ≥10x Coverage']:.2f}".rjust(10) +
        f"{stats['Heterozygous Sites']}".rjust(12) +
        f"{stats['Heterozygosity (%)']:.2f}".rjust(10) +
        f"{stats['Missing Sites']}".rjust(10) +
        f"{stats.get('Status', 'Completed')}".rjust(20)
    )
    report_lines.append(line)

report_text = "\n".join(report_lines)
with open(qc_report, "w") as f:
    f.write(report_text)

logger.info(f"✅ QC report saved to: {qc_report}")
logger.info("\nQC Report Summary:\n")
logger.info(report_text)

####################################################
# Step 6: Coverage Distribution Plots (Raw + Truncated)
####################################################
logger.info("\n===========================")
logger.info("📈 Generating coverage plots")
logger.info("===========================\n")

plot_dir = os.path.join(QC_DIR, "Plots")
os.makedirs(plot_dir, exist_ok=True)
logger.debug(f"Plot directory created: {plot_dir}")

for sample in samples:
    depth_path = depth_files[sample]
    if not os.path.exists(depth_path):
        logger.warning(f"⏭️  Skipping plot for {sample} — depth file missing.")
        continue

    logger.info(f"📊 Plotting coverage histogram for {sample}...")
    try:
        # Read raw depths
        depths_df = pd.read_csv(depth_path, sep="\t", header=None, usecols=[2], names=["depth"])
        raw_depth_values = depths_df["depth"].values

        # (A) Plot raw coverage distribution (clipped at 99th percentile for display)
        clip_thresh_raw = pd.Series(raw_depth_values).quantile(0.99)
        clipped_raw = pd.Series(raw_depth_values).clip(upper=clip_thresh_raw)

        plt.figure(figsize=(10, 5))
        plt.hist(clipped_raw, bins=100, color='skyblue', edgecolor='black')
        plt.title(f"Raw Coverage Distribution for {sample}")
        plt.xlabel(f"Depth (clipped at 99th percentile: {clip_thresh_raw:.0f})")
        plt.ylabel("Number of Positions")
        plt.grid(True)
        plt.tight_layout()

        raw_plot_path = os.path.join(plot_dir, f"{sample}_coverage_hist_raw.png")
        plt.savefig(raw_plot_path)
        plt.close()
        logger.info(f"✅ Saved raw coverage plot: {raw_plot_path}")

        # (B) Plot truncated coverage distribution
        cap = 250 if samples[sample]['data_type'] == 'ILLUMINA' else 500
        truncated_values = [min(d, cap) for d in raw_depth_values]
        # Optionally clip at 99th percentile for a nicer histogram
        clip_thresh_trunc = pd.Series(truncated_values).quantile(0.99)
        clipped_trunc = pd.Series(truncated_values).clip(upper=clip_thresh_trunc)

        plt.figure(figsize=(10, 5))
        plt.hist(clipped_trunc, bins=100, color='orange', edgecolor='black')
        plt.title(f"Truncated Coverage Distribution for {sample} (cap={cap}x)")
        plt.xlabel(f"Depth (clipped at 99th percentile: {clip_thresh_trunc:.0f})")
        plt.ylabel("Number of Positions")
        plt.grid(True)
        plt.tight_layout()

        trunc_plot_path = os.path.join(plot_dir, f"{sample}_coverage_hist_truncated.png")
        plt.savefig(trunc_plot_path)
        plt.close()
        logger.info(f"✅ Saved truncated coverage plot: {trunc_plot_path}")

    except Exception as e:
        logger.error(f"❌ Error generating coverage plots for {sample}: {e}")

logger.info("\n=======================")
logger.info("🎉 Pipeline completed!")
logger.info("=======================")
logger.info(f"Log file: {log_file_path}")
logger.info(f"QC Report: {qc_report}")
logger.info(f"Plots directory: {plot_dir}")
