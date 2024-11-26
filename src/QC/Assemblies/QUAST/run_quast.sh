#!/bin/bash

# Define paths
ASSEMBLY_DIR="/data/proj2/home/students/m.borgmann/Master_thesis/data/raw/Reference_data/assemblies"
OUTPUT_DIR="/data/proj2/home/students/m.borgmann/Master_thesis/data/raw/Quality_control/Assembly_stats/QUAST_analysis"

# First run - no minimum contig size to see distribution
mkdir -p "$OUTPUT_DIR/full_analysis"
echo "Running full analysis without contig size filter..."
quast.py \
    $ASSEMBLY_DIR/ASM2875187v1_Rapana_venosa/GCA_028751875.1_ASM2875187v1_genomic.fna \
    $ASSEMBLY_DIR/MNHN-Sthae-1_Stramonita_haemastoma/GCA_030674155.1_MNHN-Sthae-1_genomic.fna \
    $ASSEMBLY_DIR/ASM3478001v1_Urosalpinx_cinerea/GCA_034780015.1_ASM3478001v1_genomic.fna \
    $ASSEMBLY_DIR/ASM3478023v1_Concholepas_concholepas/GCA_034780235.1_ASM3478023v1_genomic.fna \
    $ASSEMBLY_DIR/ASM4096797v1_Bolinus_brandaris/GCA_040967975.1_ASM4096797v1_genomic.fna \
    $ASSEMBLY_DIR/ASM4096801v1_Hexaplex_trunculus/GCA_040968015.1_ASM4096801v1_genomic.fna \
    --threads 30 \
    --labels Rapana,Stramonita,Urosalpinx,Concholepas,Bolinus,Hexaplex \
    --eukaryote \
    --min-contig 1000 \
    --large \
    -o "$OUTPUT_DIR/full_analysis"

# Wait for first analysis to complete
echo "Analyzing contig size distribution..."
# After looking at the distribution in cumulative length plots,
# we can run the filtered analysis with a justified cutoff

# Then run with chosen cutoff based on data
# mkdir -p "$OUTPUT_DIR/filtered_analysis"
# Add second QUAST command here after reviewing first results