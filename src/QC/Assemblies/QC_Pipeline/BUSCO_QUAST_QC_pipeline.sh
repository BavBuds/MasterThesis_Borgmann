#!/bin/bash

# Define path for BUSCO lineage
BUSCO_LINEAGE="/data/proj2/home/students/m.borgmann/Master_thesis/data/raw/Quality_control/Assembly_stats/BUSCO_analysis/lineage_datasets/mollusca_odb10"

echo "Starting Assembly QC Pipeline..."
echo "================================"

# Ask for assembly file path
echo "Please enter the full path to your assembly file"
echo "Example: /data/proj2/home/students/m.borgmann/Master_thesis/data/processed/assemblies/Bolinus_SOAP/Bolinus_soap.gapClosed.fa"
echo "Path to assembly file:"
read ASSEMBLY_FILE

# Validate file exists
if [ ! -f "$ASSEMBLY_FILE" ]; then
    echo "Error: File not found: $ASSEMBLY_FILE"
    exit 1
fi

# Extract directory and name information
assembly_dir=$(dirname "$ASSEMBLY_FILE")
selected_name=$(basename "$assembly_dir")

# Create QC directories in the assembly folder
BUSCO_OUTPUT="$assembly_dir/QC_stats/${selected_name}_BUSCO"
QUAST_OUTPUT="$assembly_dir/QC_stats/${selected_name}_QUAST"

mkdir -p "$BUSCO_OUTPUT"
mkdir -p "$QUAST_OUTPUT"

echo "Using assembly: $ASSEMBLY_FILE"
echo "Output will be labeled as: $selected_name"
echo "================================"

# Run BUSCO analysis
echo "Starting BUSCO analysis..."
echo "========================="
busco -i "$ASSEMBLY_FILE" \
      -l "$BUSCO_LINEAGE" \
      -o "$selected_name" \
      -m genome \
      --out_path "$BUSCO_OUTPUT" \
      --cpu 55

# Run QUAST analysis
echo "Starting QUAST analysis..."
echo "========================="
quast.py "$ASSEMBLY_FILE" \
    --threads 55 \
    --labels "$selected_name" \
    --eukaryote \
    --min-contig 1000 \
    --large \
    -o "$QUAST_OUTPUT"

echo "QC Pipeline Completed!"
echo "===================="
echo "BUSCO results are in: $BUSCO_OUTPUT/$selected_name"
echo "QUAST results are in: $QUAST_OUTPUT"