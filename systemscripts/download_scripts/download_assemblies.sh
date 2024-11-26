#!/bin/bash

# Create a function to download and process each assembly
download_assembly() {
   local gca=$1          # GCA number needed for download path
   local assembly_id=$2   # Assembly identifier from NCBI (e.g., MNHN-Sthae-1)
   local species=$3      # Scientific name
   
   echo "Downloading $assembly_id ($species)..."
   
   # Create directory with assembly ID and species name
   mkdir -p "${assembly_id}_${species}"
   cd "${assembly_id}_${species}"
   
   wget "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/${gca:4:3}/${gca:7:3}/${gca:10:3}/${gca}_${assembly_id}/${gca}_${assembly_id}_genomic.fna.gz"
   
   # Decompress
   gunzip *_genomic.fna.gz
   
   cd ..
   echo "Completed ${assembly_id}_${species}"
}

# Make sure we're in the assemblies directory
cd /data/proj2/home/students/m.borgmann/Master_thesis/data/raw/Reference_data/assemblies

# Download each assembly using format: GCA number, Assembly ID, Species name
download_assembly "GCA_028751875.1" "ASM2875187v1" "Rapana_venosa"
download_assembly "GCA_030674155.1" "MNHN-Sthae-1" "Stramonita_haemastoma"
download_assembly "GCA_034780235.1" "ASM3478023v1" "Concholepas_concholepas"
download_assembly "GCA_040968015.1" "ASM4096801v1" "Hexaplex_trunculus"
download_assembly "GCA_034780015.1" "ASM3478001v1" "Urosalpinx_cinerea"