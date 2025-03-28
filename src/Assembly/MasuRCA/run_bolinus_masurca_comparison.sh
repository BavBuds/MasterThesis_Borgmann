#!/bin/bash

# Initialize conda for environment switching
source /data/proj2/home/students/m.borgmann/miniforge3/etc/profile.d/conda.sh

# Activate the assembly environment initially
conda activate master_thesis

# Base directory setup and configuration
BASE_DIR="/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/assemblies/Bolinus_MaSuRCA"
INPUT_DIR="/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/QCed_IlluminaFiles/Bolinus"
BUSCO_LINEAGE="/data/proj2/home/students/m.borgmann/Master_thesis/data/raw/Quality_control/Assembly_stats/BUSCO_analysis/lineage_datasets/mollusca_odb10"

# Input files
READ1="${INPUT_DIR}/SRR28863561_1.qcprocessed.fastq"
READ2="${INPUT_DIR}/SRR28863561_2.qcprocessed.fastq"

# Memory check function
check_memory() {
    local total_mem=$(free -g | awk '/^Mem:/{print $2}')
    local avail_mem=$(free -g | awk '/^Mem:/{print $7}')
    
    if [ ${avail_mem} -lt 880 ]; then  # We need at least 880GB available
        echo "ERROR: Insufficient available memory (${avail_mem}GB). Need at least 880GB"
        exit 1
    fi
    echo "Memory check passed: ${avail_mem}GB available out of ${total_mem}GB total"
}

# Logging setup
setup_logging() {
    local config_type=$1
    local timestamp=$(date +%Y%m%d_%H%M%S)
    local log_dir="${BASE_DIR}/${config_type}/logs"
    mkdir -p "${log_dir}"
    exec 3>&1 4>&2
    exec 1> >(tee -a "${log_dir}/assembly_${timestamp}.log")
    exec 2> >(tee -a "${log_dir}/assembly_${timestamp}.err" >&2)
}

# Resource monitoring function
monitor_resources() {
    local outdir=$1
    local pid=$2
    local log_file="${outdir}/logs/resource_usage.log"
    
    echo "Starting resource monitoring for PID ${pid}"
    while kill -0 $pid 2>/dev/null; do
        echo "=== Resource Usage at $(date) ===" >> "${log_file}"
        free -h >> "${log_file}"
        ps -o pid,ppid,%cpu,%mem,cmd -p $pid >> "${log_file}"
        sleep 300
    done
}

# Function to create MaSuRCA configuration
create_masurca_config() {
    local type=$1
    local outdir=$2
    
    echo "Creating MaSuRCA configuration for ${type} assembly..."
    
    # Set parameters based on assembly type
    if [ "${type}" == "permissive" ]; then
        local kmer_size=31
        local kmer_count=2
        local cgw_error=0.15
        local ovl_mer=30
        local do_trim=0
    else
        local kmer_size=61
        local kmer_count=3
        local cgw_error=0.10
        local ovl_mer=45
        local do_trim=1
    fi
    
    # Increased JF_SIZE to avoid warning
    local jf_size=1300000000

    cat > "${outdir}/masurca_config.txt" << EOF
# MaSuRCA configuration for Bolinus assembly (${type} settings)
DATA
PE= pe 150 15 ${READ1} ${READ2}
END

PARAMETERS
GRAPH_KMER_SIZE = ${kmer_size}
USE_LINKING_MATES = 1
LIMIT_JUMP_COVERAGE = 300
CA_PARAMETERS = ovlMerSize=${ovl_mer} cgwErrorRate=${cgw_error} ovlHashBits=25 ovlHashBlockLength=100000000 utgMemory=650GB obtMemory=650GB
KMER_COUNT_THRESHOLD = ${kmer_count}
NUM_THREADS = 62
JF_SIZE = ${jf_size}
DO_HOMOPOLYMER_TRIM = ${do_trim}
CLOSE_GAPS = 1
SOAP_ASSEMBLY = 0
END
EOF
}

# Function to run BUSCO and QUAST analysis
# NOTE: This function assumes we are in the QC environment.
run_qc_analysis() {
    local assembly_file=$1
    local config_type=$2
    local timestamp=$(date +%Y%m%d_%H%M%S)
    
    echo "Starting QC analysis for ${config_type} assembly at $(date)"
    
    # Run BUSCO
    echo "Running BUSCO analysis..."
    busco -i "${assembly_file}" \
          -l "${BUSCO_LINEAGE}" \
          -o "Bolinus_${config_type}_${timestamp}" \
          -m genome \
          --out_path "${BASE_DIR}/QC_stats/${config_type}/BUSCO" \
          --cpu 62
    
    # Run QUAST
    echo "Running QUAST analysis..."
    quast.py "${assembly_file}" \
             --threads 62 \
             --labels "Bolinus_${config_type}" \
             --eukaryote \
             --min-contig 10 \
             --large \
             -o "${BASE_DIR}/QC_stats/${config_type}/QUAST"
}

# Function to run MaSuRCA assembly
run_masurca() {
    local type=$1
    local outdir="${BASE_DIR}/${type}"
    local final_assembly_path="${outdir}/CA/9-terminator/genome.scf.fasta"
    local utg_file="${outdir}/CA/9-terminator/genome.utg.fasta"
    
    echo "Starting ${type} assembly at $(date)"
    mkdir -p "${outdir}"
    
    # Check if final assembly already exists
    if [ -f "${final_assembly_path}" ] && [ -s "${final_assembly_path}" ]; then
        echo "Final assembly for ${type} already exists at ${final_assembly_path}, skipping assembly."
    else
        setup_logging "${type}"
        
        # Create configuration
        create_masurca_config "${type}" "${outdir}"
        
        # Generate assembly script
        cd "${outdir}"
        echo "Generating MaSuRCA assembly script..."
        masurca masurca_config.txt
        
        if [ ! -f "assemble.sh" ]; then
            echo "ERROR: Failed to generate assembly script for ${type} configuration"
            return 1
        fi
        
        # Start assembly with resource monitoring
        echo "Starting assembly execution..."
        chmod +x assemble.sh
        ./assemble.sh &
        local assembly_pid=$!
        monitor_resources "${outdir}" ${assembly_pid} &
        local monitor_pid=$!
        
        # Wait for assembly to complete
        wait ${assembly_pid}
        local assembly_status=$?
        kill ${monitor_pid} 2>/dev/null
        
        if [ ${assembly_status} -eq 0 ] && [ -f "${final_assembly_path}" ]; then
            echo "${type} assembly completed successfully at $(date)"
        else
            echo "ERROR: ${type} assembly failed at $(date)"
            return 1
        fi
    fi
    
    # At this point, assembly is done or previously existed, we attempt QC.
    # Check if scf file is empty
    if [ ! -s "${final_assembly_path}" ]; then
        echo "No scaffolds found for ${type}. Checking unitigs..."
        if [ -s "${utg_file}" ]; then
            echo "Using unitigs at ${utg_file} for QC."
            # Switch to QC environment for QC steps
            conda deactivate
            conda activate QC
            run_qc_analysis "${utg_file}" "${type}"
            # Switch back to master_thesis environment
            conda deactivate
            conda activate master_thesis
        else
            echo "ERROR: No suitable assembly file (scaffolds or unitigs) found for QC in ${type}."
            return 1
        fi
    else
        # We have a non-empty genome.scf.fasta
        echo "Using scaffolds at ${final_assembly_path} for QC."
        # Switch to QC environment
        conda deactivate
        conda activate QC
        run_qc_analysis "${final_assembly_path}" "${type}"
        # Switch back to master_thesis environment
        conda deactivate
        conda activate master_thesis
    fi
    return 0
}

# Main execution
echo "Starting MaSuRCA assembly pipeline at $(date)"

# Check available memory
check_memory

# Create directory structure
mkdir -p "${BASE_DIR}/QC_stats"/{permissive,stringent}/{BUSCO,QUAST}
mkdir -p "${BASE_DIR}/QC_stats/comparison"

# Run assemblies sequentially

# Run permissive assembly
run_masurca "permissive"
permissive_status=$?

if [ ${permissive_status} -eq 0 ]; then
    # Run stringent assembly if permissive succeeded
    run_masurca "stringent"
    stringent_status=$?
else
    echo "Permissive assembly failed, stopping pipeline"
    exit 1
fi

# If both succeeded, attempt comparison on scf files (if they exist)
if [ ${permissive_status} -eq 0 ] && [ ${stringent_status} -eq 0 ]; then
    permissive_assembly="${BASE_DIR}/permissive/CA/9-terminator/genome.scf.fasta"
    stringent_assembly="${BASE_DIR}/stringent/CA/9-terminator/genome.scf.fasta"
    
    if [ -s "${permissive_assembly}" ] && [ -s "${stringent_assembly}" ]; then
        echo "Generating assembly comparison"
        # Switch to QC environment for QUAST comparison
        conda deactivate
        conda activate QC
        quast.py "${permissive_assembly}" \
                 "${stringent_assembly}" \
                 --threads 62 \
                 --labels "permissive,stringent" \
                 --eukaryote \
                 --large \
                 -o "${BASE_DIR}/QC_stats/comparison"
        # Switch back to master_thesis after comparison
        conda deactivate
        conda activate master_thesis
        echo "Pipeline completed successfully at $(date)"
    else
        echo "One or both assemblies lack scaffolds, cannot run final comparison."
    fi
else
    echo "Pipeline completed with errors at $(date)"
    exit 1
fi
