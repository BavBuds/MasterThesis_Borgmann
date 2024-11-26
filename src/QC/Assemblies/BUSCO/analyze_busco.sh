#!/bin/bash

# Define directories
RESULTS_DIR="/data/proj2/home/students/m.borgmann/Master_thesis/data/raw/Quality_control/Assembly_stats/BUSCO_analysis/results"
OUTPUT_DIR="/data/proj2/home/students/m.borgmann/Master_thesis/data/raw/Quality_control/Assembly_stats/BUSCO_analysis/results"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR
OUTPUT_FILE="$OUTPUT_DIR/busco_comparative_analysis.tsv"

# Function to extract metrics
extract_metrics() {
    local summary_file=$1
    local species_name=$(dirname $summary_file | xargs basename)
    
    # Extract metrics using more specific patterns
    local complete=$(grep "C:" $summary_file | sed 's/.*C:\([0-9.]*\).*/\1/')
    local single=$(grep "S:" $summary_file | sed 's/.*S:\([0-9.]*\).*/\1/')
    local duplicated=$(grep "D:" $summary_file | sed 's/.*D:\([0-9.]*\).*/\1/')
    local fragmented=$(grep "F:" $summary_file | sed 's/.*F:\([0-9.]*\).*/\1/')
    local missing=$(grep "M:" $summary_file | sed 's/.*M:\([0-9.]*\).*/\1/')
    
    printf "%-45s\t%7.1f\t%7.1f\t%7.1f\t%7.1f\t%7.1f\n" \
        "$species_name" "$complete" "$single" "$duplicated" "$fragmented" "$missing"
}

{
    # Header
    printf "%-45s\t%7s\t%7s\t%7s\t%7s\t%7s\n" \
        "Species" "Compl%" "Single%" "Dupl%" "Frag%" "Miss%"
    echo "------------------------------------------------------------------------------------------------"

    # Process files and sort by completeness
    for summary in $RESULTS_DIR/*/short_summary.specific.mollusca_odb10.*.txt; do
        extract_metrics $summary
    done | sort -k2 -nr
} | tee $OUTPUT_FILE

echo -e "\nAnalysis saved to: $OUTPUT_FILE"