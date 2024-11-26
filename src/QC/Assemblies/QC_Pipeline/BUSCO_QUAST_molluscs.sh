#!/bin/bash

# Define paths for all analyses
ASSEMBLIES_DIR="/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/assemblies"
BUSCO_LINEAGE="/data/proj2/home/students/m.borgmann/Master_thesis/data/raw/Quality_control/Assembly_stats/BUSCO_analysis/lineage_datasets/mollusca_odb10"
QC_OUTPUT_BASE="/data/proj2/home/students/m.borgmann/Master_thesis/data/raw/Quality_control/Assembly_stats"

# Create necessary output directories
BUSCO_OUTPUT="${QC_OUTPUT_BASE}/BUSCO_analysis/results"
QUAST_OUTPUT="${QC_OUTPUT_BASE}/QUAST_analysis"
mkdir -p "$BUSCO_OUTPUT"
mkdir -p "$QUAST_OUTPUT"

echo "Starting Assembly QC Pipeline..."
echo "================================"

# List available assemblies and let user choose
echo "Available assemblies:"
echo "-------------------"
# Get directories and store in array
available_assemblies=($(ls -d ${ASSEMBLIES_DIR}/*))
assembly_names=($(ls -d ${ASSEMBLIES_DIR}/* | xargs -n 1 basename))

# Print numbered list of assemblies
for i in "${!assembly_names[@]}"; do
    echo "$((i+1))) ${assembly_names[$i]}"
done

# Get user selection
echo ""
echo "Enter the number of the assembly you want to analyze:"
read selection

# Validate input
if ! [[ "$selection" =~ ^[0-9]+$ ]] || [ "$selection" -lt 1 ] || [ "$selection" -gt "${#assembly_names[@]}" ]; then
    echo "Invalid selection. Please run script again and select a number between 1 and ${#assembly_names[@]}."
    exit 1
fi

# Get selected assembly
selected_idx=$((selection-1))
selected_assembly="${available_assemblies[$selected_idx]}"
selected_name="${assembly_names[$selected_idx]}"

echo "You selected: $selected_name"
echo "================================"

# Function to find the assembly file
find_assembly_file() {
    local assembly_dir=$1
    # Look for common assembly file extensions
    for ext in ".gapclosed.fasta" ".fasta" ".fa" ".scaf" ".scaffold.fasta" ".final.fasta"; do
        local file=$(find "$assembly_dir" -maxdepth 1 -name "*$ext" | head -n 1)
        if [ ! -z "$file" ]; then
            echo "$file"
            return
        fi
    done
    echo ""
}

# Find the assembly file
ASSEMBLY_FILE=$(find_assembly_file "$selected_assembly")

if [ -z "$ASSEMBLY_FILE" ]; then
    echo "Error: Could not find assembly file in $selected_assembly"
    exit 1
fi

echo "Found assembly file: $ASSEMBLY_FILE"
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
    -o "$QUAST_OUTPUT/${selected_name}_analysis"

echo "QC Pipeline Completed!"
echo "===================="
echo "BUSCO results are in: $BUSCO_OUTPUT/$selected_name"
echo "QUAST results are in: $QUAST_OUTPUT/${selected_name}_analysis"