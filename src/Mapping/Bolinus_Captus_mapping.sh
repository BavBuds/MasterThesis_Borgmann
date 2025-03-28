#!/bin/bash

# Define input and output paths
GENOME="/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/assemblies/Bolinus_Captus/SRR28863561__captus-asm/02_assemblies/assembly.fasta"
READ1="/data/proj2/home/students/m.borgmann/Master_thesis/data/raw/Reference_data/WGS_projects/PRJNA1106534_Bolinus_brandaris/SRR28863561/SRR28863561_R1.fastq"
READ2="/data/proj2/home/students/m.borgmann/Master_thesis/data/raw/Reference_data/WGS_projects/PRJNA1106534_Bolinus_brandaris/SRR28863561/SRR28863561_R2.fastq"
MAPPING_DIR="/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/assemblies/Bolinus_Captus/Mapping"

# Create output directory if it doesn't exist
echo "Creating output directory: $MAPPING_DIR"
mkdir -p "$MAPPING_DIR"

# Index the genome
echo "Indexing the genome assembly..."
bwa index "$GENOME"

# Run BWA-MEM alignment
echo "Running BWA-MEM alignment..."
bwa mem -t 35 \ 
  "$GENOME" \ 
  "$READ1" \ 
  "$READ2" \ 
  > "$MAPPING_DIR/mapped_reads.sam"

# Convert SAM to BAM and sort
echo "Converting SAM to BAM and sorting..."
samtools view -Sb "$MAPPING_DIR/mapped_reads.sam" | \
  samtools sort -o "$MAPPING_DIR/mapped_reads.sorted.bam"

# Index the BAM file
echo "Indexing the BAM file..."
samtools index "$MAPPING_DIR/mapped_reads.sorted.bam"

# Cleanup: remove intermediate SAM file
echo "Cleaning up intermediate files..."
rm "$MAPPING_DIR/mapped_reads.sam"

echo "Mapping pipeline completed successfully. Outputs saved in $MAPPING_DIR"
