#!/usr/bin/env bash
###############################################################################
# Script: create_20mb_subset_multithreaded.sh
# Description:
#   1) Sort contigs by length (descending) using multithreading.
#   2) Select the largest contigs until reaching ~20 Mb total length.
#   3) Extract those contigs into a new FASTA (subset_20mb.fasta) using multithreading.
#   4) Run a quality check on the subset using parallel coverage analysis.
#   5) Output everything into the PSMC directory.
###############################################################################

set -euo pipefail

# -----------------------------
# 0) User-defined paths
# -----------------------------
ASSEMBLY="/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/assemblies/Bolinus_Captus/SRR28863561__captus-asm/02_assemblies/assembly.fasta"
BAM="/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/assemblies/Bolinus_Captus/Mapping/mapped_reads.sorted.bam"

# Output directory
PSMC_DIR="/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/Popgen_analysis/PSMC"
mkdir -p "${PSMC_DIR}"

# Target genome size for subset (in bp) - 20,000,000 for 20 Mb
TARGET_BP=20000000

# File paths we'll create
SORTED_CONTIGS="${PSMC_DIR}/sorted_contigs.txt"
CONTIGS_20MB_TXT="${PSMC_DIR}/contigs_20mb.txt"
SUBSET_FA="${PSMC_DIR}/subset_20mb.fasta"

# Number of threads
THREADS=40

echo "Starting creation of 20 Mb subset..."
date
echo ""

# -----------------------------
# 1) Sort contigs by length (descending)
# -----------------------------
echo "Sorting contigs by length..."
sort --parallel="${THREADS}" -k2,2nr "${ASSEMBLY}.fai" > "${SORTED_CONTIGS}"
echo "Sorted contigs written to: ${SORTED_CONTIGS}"

# -----------------------------
# 2) Select contigs until we reach ~20 Mb
# -----------------------------
echo "Selecting contigs until we reach ~20 Mb..."
awk -v target="${TARGET_BP}" '
BEGIN { sum=0 }
{
  if (sum < target) {
    print $1;  # Print only the contig name
    sum += $2;  # Add the contig length to the running total
  }
}
' "${SORTED_CONTIGS}" > "${CONTIGS_20MB_TXT}"
echo "Contigs used for the 20 Mb subset written to: ${CONTIGS_20MB_TXT}"

# -----------------------------
# 3) Extract those contigs to a new FASTA
# -----------------------------
echo "Extracting selected contigs to FASTA using multithreading..."
mkdir -p temp_fasta_files
cat "${CONTIGS_20MB_TXT}" | xargs -P "${THREADS}" -I {} sh -c \
"samtools faidx '${ASSEMBLY}' {} > temp_fasta_files/{}.fa"

# Combine the individual FASTA files
cat temp_fasta_files/*.fa > "${SUBSET_FA}"
rm -r temp_fasta_files
echo "20 Mb subset FASTA created at: ${SUBSET_FA}"

# -----------------------------
# 4) Expanded Quality Check
# -----------------------------
echo ""
echo "=== Expanded Quality Analysis ==="

# Calculate total length of the subset FASTA
echo "Calculating total length of the subset FASTA..."
SUBSET_LEN=$(awk '/^>/ {next} {sum += length($0)} END {print sum}' "${SUBSET_FA}")
echo "  - Subset length: ${SUBSET_LEN} bp"

# Count missing data (Ns)
echo "Counting missing data (N's) in the subset..."
TOTAL_NS=$(grep -o "N" "${SUBSET_FA}" | wc -l || echo 0)
MISSING_PCT=$(echo "scale=2; (${TOTAL_NS} / ${SUBSET_LEN}) * 100" | bc)
echo "  - Total Ns: ${TOTAL_NS}"
echo "  - Missing data: ${MISSING_PCT}%"

# Coverage analysis: generate depth file for selected contigs
COV_FILE="${PSMC_DIR}/subset_20mb_depth.txt"
mkdir -p temp_depth_files
echo "Calculating per-base coverage for the selected contigs using multithreading..."
cat "${CONTIGS_20MB_TXT}" | xargs -P "${THREADS}" -I {} sh -c \
"samtools depth -r {} -q 0 -Q 0 '${BAM}' > temp_depth_files/{}.depth"

# Combine the individual depth files
cat temp_depth_files/*.depth > "${COV_FILE}"
rm -r temp_depth_files

# Mean coverage calculation
echo "Calculating mean coverage across all sites..."
if [[ -s "${COV_FILE}" ]]; then
    MEAN_COV=$(awk '{ sum+=$3; count++ } END { if(count>0) print sum/count; else print 0 }' "${COV_FILE}")
else
    MEAN_COV=0
fi
echo "  - Mean coverage: ${MEAN_COV}X"

# Calculate fraction of sites covered at least 10x
echo "Calculating fraction of sites with coverage >= 10x..."
if [[ -s "${COV_FILE}" ]]; then
    TOTAL_SITES=$(awk 'END {print NR}' "${COV_FILE}")
    SITES_10X=$(awk '$3 >= 10 {count++} END {print count}' "${COV_FILE}")
    FRAC_10X=$(echo "scale=2; (${SITES_10X} / ${TOTAL_SITES}) * 100" | bc)
else
    TOTAL_SITES=0
    SITES_10X=0
    FRAC_10X=0
fi
echo "  - Total sites: ${TOTAL_SITES}"
echo "  - Sites with >=10x coverage: ${SITES_10X} (${FRAC_10X}%)"

# Summary report
echo ""
echo "=== Quality Summary ==="
echo "  - Total Length: ${SUBSET_LEN} bp"
echo "  - Missing Data: ${MISSING_PCT}%"
echo "  - Mean Coverage: ${MEAN_COV}X"
echo "  - Sites with >=10x Coverage: ${SITES_10X} (${FRAC_10X}%)"
echo ""
