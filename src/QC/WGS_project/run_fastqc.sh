#!/bin/bash

# Define paths
DATA_DIR="/data/proj2/home/students/m.borgmann/Master_thesis/data/raw/Reference_data/WGS_projects"
OUTPUT_DIR="/data/proj2/home/students/m.borgmann/Master_thesis/data/raw/Quality_control/Fastq_stats"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Function to run FastQC on paired reads
run_fastqc() {
    local species_dir=$1
    local srr_number=$2
    
    echo "Processing $species_dir - $srr_number"
    
    # Create species-specific output directory
    local species_output="$OUTPUT_DIR/$(basename $species_dir)"
    mkdir -p "$species_output"
    
    # Run FastQC with multiple threads
    fastqc \
        "$DATA_DIR/$species_dir/$srr_number/${srr_number}_1.fastq" \
        "$DATA_DIR/$species_dir/$srr_number/${srr_number}_2.fastq" \
        --outdir="$species_output" \
        --threads 30
}

# Run FastQC for each dataset
run_fastqc "PRJNA1106534_Bolinus_brandaris" "SRR28863561"
run_fastqc "PRJNA1106542_Hexaplex_trunculus" "SRR28865916"

echo "FastQC analysis completed. Results are in: $OUTPUT_DIR"