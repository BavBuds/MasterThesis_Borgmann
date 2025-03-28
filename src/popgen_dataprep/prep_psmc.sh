#!/usr/bin/env bash
###############################################################################
# Script: prep_psmc.sh
# Purpose:
#   1) Create a PSMC input file from a subset of contigs (using >10x coverage).
#   2) Run PSMC on that dataset.
#   3) Perform 100 bootstrap replicates for confidence intervals.
#
# Author: Your Name
# Date: YYYY-MM-DD
###############################################################################

set -euo pipefail

#######################
# 0) User-defined paths
#######################

# Input files
REFERENCE="/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/assemblies/Bolinus_Captus/SRR28863561__captus-asm/02_assemblies/assembly.fasta"
BAM="/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/assemblies/Bolinus_Captus/Mapping/mapped_reads.sorted.bam"
CONTIGS_20MB_TXT="/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/Popgen_analysis/PSMC/dataprep/contigs_20mb.txt"

# Intermediate files directory
DATAPREP_DIR="/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/Popgen_analysis/PSMC/dataprep"
mkdir -p "${DATAPREP_DIR}"

# Final PSMC outputs directory
PSMC_DIR="/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/Popgen_analysis/PSMC"
mkdir -p "${PSMC_DIR}"

# Number of threads (adjust as needed)
THREADS=62

# Coverage thresholds for vcf2fq
MIN_COV=10
MAX_COV=999999  # you can set a max if you want to exclude extremely high coverage

#############################
# 1) Generate diploid FASTQ
#############################

echo ">>> Generating diploid FASTQ with coverage >= ${MIN_COV}..."

# We’ll use samtools/bcftools to:
#   (a) mpileup restricted to subset contigs
#   (b) call variants
#   (c) convert to FASTQ format with coverage filters
# If you used bcftools to do advanced filtering before, you can adapt that part here.
# Note that -B disables BAQ computation (for PSMC, recommended).
# -Q sets minimum base quality for mpileup to consider.

bcftools mpileup \
    --threads "${THREADS}" \
    -Q 20 \
    -f "${REFERENCE}" \
    -R "${CONTIGS_20MB_TXT}" \
    -B \
    "${BAM}" \
| bcftools call \
    --threads "${THREADS}" \
    -c \
| vcfutils.pl vcf2fq \
    -d "${MIN_COV}" \
    -D "${MAX_COV}" \
> "${DATAPREP_DIR}/subset_20mb.diploid.fq"

echo ">>> Diploid FASTQ created: ${DATAPREP_DIR}/subset_20mb.diploid.fq"

###################################
# 2) Convert FASTQ to PSMC input
###################################

echo ">>> Converting FASTQ to PSMCFA..."

fq2psmcfa \
    -q20 \
    "${DATAPREP_DIR}/subset_20mb.diploid.fq" \
    > "${DATAPREP_DIR}/subset_20mb.psmcfa"

echo ">>> PSMCFA created: ${DATAPREP_DIR}/subset_20mb.psmcfa"

#################################
# 3) Run PSMC on real dataset
#################################

echo ">>> Running PSMC on real dataset..."

# You can adjust the PSMC parameters as needed:
#   -N25 : max 25 iterations
#   -t15 : initial theta
#   -r5  : initial rho
#   -p 4+25*2+4+6 : recommended pattern
psmc \
    -N25 \
    -t15 \
    -r5 \
    -p "4+25*2+4+6" \
    -o "${PSMC_DIR}/subset_20mb.psmc" \
    "${DATAPREP_DIR}/subset_20mb.psmcfa"

echo ">>> PSMC run complete: ${PSMC_DIR}/subset_20mb.psmc"

#####################################################
# 4) Perform 100 bootstrap replicates for confidence
#####################################################

echo ">>> Starting 100 bootstrap replicates..."

# We'll create a subfolder for bootstrap results
BOOTSTRAP_DIR="${PSMC_DIR}/bootstrap"
mkdir -p "${BOOTSTRAP_DIR}"

# Typical approach:
#  1) For each replicate i in 1..100:
#     a) "splitfa" randomly re-samples blocks from the .psmcfa
#     b) run psmc on that re-sampled data
#     c) output each run to a separate file for later analysis

for i in $(seq 1 100); do
    BOOT_PSMCFA="${BOOTSTRAP_DIR}/bootstrap_${i}.psmcfa"
    BOOT_OUT="${BOOTSTRAP_DIR}/bootstrap_${i}.psmc"

    echo "  - Bootstrap replicate ${i}..."

    # Re-sample segments from the original PSMCFA
    splitfa "${DATAPREP_DIR}/subset_20mb.psmcfa" > "${BOOT_PSMCFA}"

    # Run PSMC with the same parameters as the real dataset
    psmc \
        -N25 \
        -t15 \
        -r5 \
        -p "4+25*2+4+6" \
        -o "${BOOT_OUT}" \
        "${BOOT_PSMCFA}"
done

echo ">>> Bootstrapping complete. Results in: ${BOOTSTRAP_DIR}"

############################################
# 5) (Optional) Summarize or plot results
############################################
# Common practice is to plot the original PSMC curve along with
#  the 100 bootstrap curves for confidence intervals, e.g. with
#  the 'psmc_plot.pl' script or custom R/Python code.
# 
# For example:
#   psmc_plot.pl -p "Bolinus_captus" Bolinus_captus_psmc_plot subset_20mb.psmc \
#                bootstrap/bootstrap_*.psmc
#
# You can do that manually or add here as needed.

echo ""
echo "=== All steps complete. ==="
echo "Final PSMC output is in:    ${PSMC_DIR}/subset_20mb.psmc"
echo "Bootstrap files are in:     ${BOOTSTRAP_DIR}"
echo "PSMC input .psmcfa is here: ${DATAPREP_DIR}/subset_20mb.psmcfa"
