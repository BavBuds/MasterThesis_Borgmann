#!/bin/bash

# Check if assembly file path and kmer suffix are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <assembly_file_path> <kmer_suffix>"
    exit 1
fi

ASSEMBLY_FILE=$1
KMER_SUFFIX=$2

# Define paths
BUSCO_LINEAGE="/data/proj2/home/students/m.borgmann/Master_thesis/data/raw/Quality_control/Assembly_stats/BUSCO_analysis/lineage_datasets/mollusca_odb10"
BASE_DIR="/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/assemblies/Bolinus_SOAP"

# Create output directories
BUSCO_OUTPUT="$BASE_DIR/QC_stats/${KMER_SUFFIX}_BUSCO"
QUAST_OUTPUT="$BASE_DIR/QC_stats/${KMER_SUFFIX}_QUAST"

mkdir -p "$BUSCO_OUTPUT"
mkdir -p "$QUAST_OUTPUT"

echo "Starting QC Pipeline for k$KMER_SUFFIX assembly..."
echo "================================"

# Run BUSCO analysis
echo "Starting BUSCO analysis..."
echo "========================="
busco -i "$ASSEMBLY_FILE" \
      -l "$BUSCO_LINEAGE" \
      -o "Bolinus_SOAP_${KMER_SUFFIX}" \
      -m genome \
      --out_path "$BUSCO_OUTPUT" \
      --cpu 55

# Run QUAST analysis
echo "Starting QUAST analysis..."
echo "========================="
quast.py "$ASSEMBLY_FILE" \
    --threads 55 \
    --labels "Bolinus_SOAP_${KMER_SUFFIX}" \
    --eukaryote \
    --min-contig 1000 \
    --large \
    -o "$QUAST_OUTPUT"

echo "QC Pipeline Completed for k$KMER_SUFFIX!"
echo "===================="