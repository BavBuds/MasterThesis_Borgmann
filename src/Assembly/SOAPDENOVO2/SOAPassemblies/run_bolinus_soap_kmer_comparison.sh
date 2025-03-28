#!/bin/bash

# Base directories
BASE_DIR="/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/assemblies"
SOAP_DIR="$BASE_DIR/Bolinus_SOAP"

# Create directory structure
mkdir -p "$SOAP_DIR/k63"
mkdir -p "$SOAP_DIR/k127"
mkdir -p "$SOAP_DIR/QC_stats/k63_BUSCO"
mkdir -p "$SOAP_DIR/QC_stats/k63_QUAST"
mkdir -p "$SOAP_DIR/QC_stats/k127_BUSCO"
mkdir -p "$SOAP_DIR/QC_stats/k127_QUAST"
mkdir -p "$SOAP_DIR/QC_stats/comparison"
mkdir -p "$SOAP_DIR/QC_stats/comparison/BUSCO"
mkdir -p "$SOAP_DIR/QC_stats/comparison/QUAST"

# Input files
READ1="/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/QCed_IlluminaFiles/Bolinus/SRR28863561_1.qcprocessed.fastq"
READ2="/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/QCed_IlluminaFiles/Bolinus/SRR28863561_2.qcprocessed.fastq"

# Function to create config and run SOAP
run_soap() {
    local kmer=$1
    local outdir="$SOAP_DIR/k${kmer}"
    local prefix="$outdir/Bolinus_soap_k${kmer}"
    
    # Create config file
    cat > "$outdir/soap_config.txt" << EOF
max_rd_len=150
[LIB]
avg_ins=150
reverse_seq=0
asm_flags=3
rank=1
pair_num_cutoff=5
map_len=40
q1=$READ1
q2=$READ2
EOF

    echo "Starting SOAPdenovo2 assembly with K=$kmer..."
    
    # Run assembly
    SOAPdenovo-${kmer}mer pregraph \
        -s "$outdir/soap_config.txt" \
        -K $kmer \
        -p 55 \
        -o $prefix \
        > "$outdir/pregraph.log" 2>&1

    SOAPdenovo-${kmer}mer contig \
        -g $prefix \
        > "$outdir/contig.log" 2>&1

    SOAPdenovo-${kmer}mer map \
        -s "$outdir/soap_config.txt" \
        -g $prefix \
        -p 55 \
        > "$outdir/map.log" 2>&1

    SOAPdenovo-${kmer}mer scaff \
        -g $prefix \
        -F \
        > "$outdir/scaff.log" 2>&1

    GapCloser \
        -b "$outdir/soap_config.txt" \
        -a "${prefix}.scafSeq" \
        -o "${prefix}.gapclosed.fasta" \
        -t 55 \
        > "$outdir/gapcloser.log" 2>&1
}

# Run both assemblies
run_soap 63
run_soap 127

# Define QC pipeline script path
QC_SCRIPT="/data/proj2/home/students/m.borgmann/Master_thesis/src/QC/Assemblies/QC_Pipeline/qc_pipeline_kmertrial_BolinusSOAPdenovo2.sh"

# Run QC pipeline for both assemblies
echo "Running QC pipeline for k63 assembly..."
$QC_SCRIPT "$SOAP_DIR/k63/Bolinus_soap_k63.gapclosed.fasta" "k63"

echo "Running QC pipeline for k127 assembly..."
$QC_SCRIPT "$SOAP_DIR/k127/Bolinus_soap_k127.gapclosed.fasta" "k127"