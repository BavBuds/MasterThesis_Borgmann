#!/bin/bash
# Create output directory for Bolinus SOAP assembly
OUTPUT_DIR="/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/assemblies/Bolinus_SOAP"
mkdir -p $OUTPUT_DIR

# Input files
READ1="/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/QCed_IlluminaFiles/Bolinus/SRR28863561_1.qcprocessed.fastq"
READ2="/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/QCed_IlluminaFiles/Bolinus/SRR28863561_2.qcprocessed.fastq"

# Create config file
CONFIG_FILE="$OUTPUT_DIR/soap_config.txt"
cat > $CONFIG_FILE << EOF
max_rd_len=150
[LIB]
avg_ins=300
reverse_seq=0
asm_flags=3
rank=1
pair_num_cutoff=3
map_len=32
q1=$READ1
q2=$READ2
EOF

# Set parameters for SOAPdenovo2
K=127 # K-mer size
N_THREADS=55 # Using 55 cores
PREFIX="$OUTPUT_DIR/Bolinus_soap"

# Run SOAPdenovo2
echo "Starting SOAPdenovo2 assembly for Bolinus with K=$K..."

# Pregraph step
SOAPdenovo-127mer pregraph \
    -s $CONFIG_FILE \
    -K $K \
    -p $N_THREADS \
    -o $PREFIX \
    > $OUTPUT_DIR/pregraph.log 2>&1

# Contig construction
SOAPdenovo-127mer contig \
    -g $PREFIX \
    > $OUTPUT_DIR/contig.log 2>&1

# Map reads to contigs
SOAPdenovo-127mer map \
    -s $CONFIG_FILE \
    -g $PREFIX \
    -p $N_THREADS \
    > $OUTPUT_DIR/map.log 2>&1

# Scaffolding
SOAPdenovo-127mer scaff \
    -g $PREFIX \
    -F \
    > $OUTPUT_DIR/scaff.log 2>&1

# Gap filling
GapCloser \
    -b $CONFIG_FILE \
    -a ${PREFIX}.scafSeq \
    -o ${PREFIX}.gapclosed.fasta \
    -t $N_THREADS \
    > $OUTPUT_DIR/gapcloser.log 2>&1

echo "Assembly completed. Check logs in $OUTPUT_DIR for any errors."
echo "Final assembly is in ${PREFIX}.gapclosed.fasta"

# Calculate basic statistics
echo "Calculating basic assembly statistics..."
grep -v "^>" ${PREFIX}.gapclosed.fasta | tr -d '\n' | wc -c > $OUTPUT_DIR/assembly_length.txt
grep "^>" ${PREFIX}.gapclosed.fasta | wc -l > $OUTPUT_DIR/scaffold_count.txt

echo "Done! Assembly and statistics are in $OUTPUT_DIR"