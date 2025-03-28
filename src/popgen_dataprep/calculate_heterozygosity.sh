#!/usr/bin/env bash
###############################################################################
# Script: calculate_heterozygosity_v1.sh
# Description:
#   1) Use bcftools mpileup on the *full* assembly, restricted to the subset 
#      contigs listed in contigs_20mb.txt.
#   2) Call SNPs with bcftools call.
#   3) Filter variants by depth, quality, etc.
#   4) Calculate heterozygosity as SNPs / total bp in the 20 Mb subset.
###############################################################################

set -euo pipefail

# -----------------------------
# 0) User-defined paths
# -----------------------------
# The full assembly used in the original mapping
FULL_ASSEMBLY="/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/assemblies/Bolinus_Captus/SRR28863561__captus-asm/02_assemblies/assembly.fasta"

# Your original BAM aligned against the full assembly
BAM="/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/assemblies/Bolinus_Captus/Mapping/mapped_reads.sorted.bam"

# Text file listing the subset contigs (one contig name per line)
CONTIGS_20MB_TXT="/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/Popgen_analysis/PSMC/dataprep/contigs_20mb.txt"

# The already-created 20 Mb subset FASTA (for length calculation)
SUBSET_FASTA="/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/Popgen_analysis/PSMC/dataprep/subset_20mb.fasta"

# Where to place VCFs and final heterozygosity results
OUTPUT_DIR="/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/Popgen_analysis/PSMC/dataprep"

# Number of threads (adjust as needed)
THREADS=60

# -----------------------------
# 1) Define output file names
# -----------------------------
MPILEUP_VCF="${OUTPUT_DIR}/subset_20mb.mpileup.vcf.gz"
RAW_CALLS_VCF="${OUTPUT_DIR}/subset_20mb.raw_calls.vcf.gz"
FILTERED_VCF="${OUTPUT_DIR}/subset_20mb.filtered.vcf.gz"
HET_RESULTS="${OUTPUT_DIR}/heterozygosity_results.txt"

# Make sure the output directory exists
mkdir -p "${OUTPUT_DIR}"

# -----------------------------
# 2) Generate mpileup (raw genotype likelihoods)
# -----------------------------
echo ">>> Generating mpileup restricted to subset contigs..."

bcftools mpileup \
    --threads "${THREADS}" \
    -f "${FULL_ASSEMBLY}" \
    -R "${CONTIGS_20MB_TXT}" \
    "${BAM}" \
    -Oz \
    -o "${MPILEUP_VCF}"

bcftools index "${MPILEUP_VCF}"
echo ">>> Mpileup complete: ${MPILEUP_VCF}"

# -----------------------------
# 3) Call variants
# -----------------------------
echo ">>> Calling variants..."
bcftools call \
    --threads "${THREADS}" \
    -mv \
    -Oz \
    -o "${RAW_CALLS_VCF}" \
    "${MPILEUP_VCF}"

bcftools index "${RAW_CALLS_VCF}"
echo ">>> Variant calling complete: ${RAW_CALLS_VCF}"

# -----------------------------
# 4) Filter variants
# -----------------------------
echo ">>> Filtering variants..."
bcftools filter \
    -i "DP>=10 && DP<=80 && QUAL>=30" \
    -Oz \
    -o "${FILTERED_VCF}" \
    "${RAW_CALLS_VCF}"

bcftools index "${FILTERED_VCF}"
echo ">>> Filtering complete: ${FILTERED_VCF}"

# -----------------------------
# 5) Calculate heterozygosity
# -----------------------------
echo ">>> Calculating heterozygosity..."

# Count total SNPs
SNP_COUNT=$(bcftools view -H "${FILTERED_VCF}" | wc -l)

# Calculate total bp in the 20 Mb subset
TOTAL_BP=$(awk '/^>/ {next} {sum += length($0)} END {print sum}' "${SUBSET_FASTA}")

# Heterozygosity = SNPs per bp
HET=$(echo "scale=6; ${SNP_COUNT} / ${TOTAL_BP}" | bc)

# Save results
{
  echo "=== Heterozygosity Results ==="
  echo "Total base pairs (subset): ${TOTAL_BP}"
  echo "Total SNPs: ${SNP_COUNT}"
  echo "Heterozygosity (SNPs / bp): ${HET}"
} > "${HET_RESULTS}"

echo ""
echo ">>> Heterozygosity calculation complete."
cat "${HET_RESULTS}"
