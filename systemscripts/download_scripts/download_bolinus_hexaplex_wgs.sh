#!/bin/bash

# Author: Your Name
# Date: 2024-10-18
# Description: Script to download SRA data for Bolinus brandaris and Hexaplex trunculus using SRA Toolkit and fastq-dump.

# Exit immediately if a command exits with a non-zero status
set -e

# Directory where SRA data will be downloaded and processed
BASE_DIR="/data/proj2/home/students/m.borgmann/Master_thesis/data/raw/reference_genomes/WGS_projects"

# Create directories for each project
mkdir -p "$BASE_DIR/PRJNA1106534"  # Bolinus brandaris
mkdir -p "$BASE_DIR/PRJNA1106601"  # Hexaplex trunculus

# Bolinus brandaris: Download SRA data
echo "Downloading SRA data for Bolinus brandaris (PRJNA1106534)..."
prefetch --max-size 50G --output-directory "$BASE_DIR/PRJNA1106534" PRJNA1106534

# Bolinus brandaris: Convert .sra files to FASTQ
echo "Converting SRA data for Bolinus brandaris to FASTQ..."
fastq-dump --split-3 --outdir "$BASE_DIR/PRJNA1106534" "$BASE_DIR/PRJNA1106534"/*.sra

# Hexaplex trunculus: Download SRA data
echo "Downloading SRA data for Hexaplex trunculus (PRJNA1106601)..."
prefetch --max-size 50G --output-directory "$BASE_DIR/PRJNA1106601" PRJNA1106601

# Hexaplex trunculus: Convert .sra files to FASTQ
echo "Converting SRA data for Hexaplex trunculus to FASTQ..."
fastq-dump --split-3 --outdir "$BASE_DIR/PRJNA1106601" "$BASE_DIR/PRJNA1106601"/*.sra

echo "Download and conversion process completed for both species."
