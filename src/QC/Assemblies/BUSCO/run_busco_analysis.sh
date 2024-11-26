#!/bin/bash

# Save as /data/proj2/home/students/m.borgmann/Master_thesis/src/run_busco_analysis.sh

# Define paths
ASSEMBLY_DIR="/data/proj2/home/students/m.borgmann/Master_thesis/data/raw/Reference_data/assemblies"
LINEAGE_PATH="/data/proj2/home/students/m.borgmann/Master_thesis/data/raw/Quality_control/Assembly_stats/BUSCO_analysis/lineage_datasets/mollusca_odb10"
OUTPUT_DIR="/data/proj2/home/students/m.borgmann/Master_thesis/data/raw/Quality_control/Assembly_stats/BUSCO_analysis/results"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Function to run BUSCO
run_busco() {
    local assembly_dir=$1
    local species_name=$(basename $assembly_dir)
    local fna_file=$(find $assembly_dir -name "*genomic.fna")
    
    echo "Processing $species_name..."
    
    busco -i $fna_file \
          -l $LINEAGE_PATH \
          -o $species_name \
          -m genome \
          --out_path $OUTPUT_DIR \
          --cpu 35
}

# Process each assembly
for assembly_dir in $ASSEMBLY_DIR/*; do
    run_busco $assembly_dir
done