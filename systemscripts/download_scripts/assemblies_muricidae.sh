#!/bin/bash

# Set the output directory
OUTPUT_DIR="/data/proj2/home/students/m.borgmann/Master_thesis/data/reference_genome/assemblies"

# Array of assembly accessions to download
ASSEMBLIES=(
    "GCA_028751875.1"
    "GCA_030674155.1"
    "GCA_034780235.1"
    "GCA_040968015.1"
    "GCA_034780015.1"
)

# Function to download and process an assembly
download_assembly() {
    accession=$1
    
    echo "Downloading assembly: $accession"
    
    # Use ncbi-datasets to download the assembly
    datasets download genome accession $accession --dehydrated --include genome --filename $OUTPUT_DIR/$accession.zip
    
    # Extract the downloaded assembly
    unzip -d $OUTPUT_DIR/$accession $OUTPUT_DIR/$accession.zip
    rm $OUTPUT_DIR/$accession.zip  # Remove the zip file after extraction
    
    echo "Finished processing $accession"
}

# Main execution
mkdir -p $OUTPUT_DIR

for assembly in "${ASSEMBLIES[@]}"; do
    download_assembly $assembly
done

echo "All assemblies have been downloaded and processed."
