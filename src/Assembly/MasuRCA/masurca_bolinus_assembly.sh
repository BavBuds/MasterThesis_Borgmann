#!/bin/bash
# Create output directory for Bolinus MaSuRCA assembly
OUTPUT_DIR="/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/assemblies/Bolinus_MaSuRCA"
mkdir -p $OUTPUT_DIR

# Input files
READ1="/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/QCed_IlluminaFiles/Bolinus/SRR28863561_1.qcprocessed.fastq"
READ2="/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/QCed_IlluminaFiles/Bolinus/SRR28863561_2.qcprocessed.fastq"

# Create MaSuRCA configuration file
CONFIG_FILE="$OUTPUT_DIR/masurca_config.txt"
cat > $CONFIG_FILE << EOF
# MaSuRCA configuration for Bolinus assembly
DATA
PE= pe 150 15 $READ1 $READ2
END

PARAMETERS
# Adjusted for ~2Gb genome
GRAPH_KMER_SIZE = 31
USE_LINKING_MATES = 1
LIMIT_JUMP_COVERAGE = 60
CA_PARAMETERS = ovlMerSize=30 cgwErrorRate=0.15 ovlMemory=4GB
KMER_COUNT_THRESHOLD = 2
NUM_THREADS = 38
# Increased for 2Gb genome size
JF_SIZE = 20000000000
DO_HOMOPOLYMER_TRIM = 0
CLOSE_GAPS = 1
SOAP_ASSEMBLY = 0
END
EOF

# Navigate to output directory
cd $OUTPUT_DIR

# Generate assembly script
masurca $CONFIG_FILE

# Start assembly
echo "Starting MaSuRCA assembly for Bolinus..."
start_time=$(date +%s)
./assemble.sh

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
echo "Assembly took $((elapsed_time/3600)) hours, $((elapsed_time%3600/60)) minutes"
echo "Assembly completed. Check final output in $OUTPUT_DIR/CA"