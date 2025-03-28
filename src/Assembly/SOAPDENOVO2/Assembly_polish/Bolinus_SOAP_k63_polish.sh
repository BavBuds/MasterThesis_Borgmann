#!/bin/bash
set -e  # Exit immediately if a command exits with a non-zero status

# Redirect all output and error messages to a log file
LOG_FILE="/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/assemblies/Bolinus_SOAP/k63_polish/pipeline.log"
exec > >(tee -i "$LOG_FILE") 2>&1

echo "Pipeline started on $(date)"

# Directories
ASSEMBLY_DIR="/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/assemblies/Bolinus_SOAP/k63"
OUTPUT_DIR="/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/assemblies/Bolinus_SOAP/k63_polish"
POLISHED_ASSEMBLY_DIR="$OUTPUT_DIR/Bolinus_k63_polished"
MAPPING_DIR="$OUTPUT_DIR/Mapping_files"
QC_DIR="/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/assemblies/Bolinus_SOAP/QC_stats"
QUAST_OUTPUT="$QC_DIR/k63_polished_QUAST"
BUSCO_OUTPUT="$QC_DIR/k63_polished_BUSCO"

# Files
ASSEMBLY_FILE="$ASSEMBLY_DIR/Bolinus_soap_k63.gapclosed.fasta"
READ1="/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/QCed_IlluminaFiles/Bolinus/SRR28863561_1.qcprocessed.fastq"
READ2="/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/QCed_IlluminaFiles/Bolinus/SRR28863561_2.qcprocessed.fastq"
POLISHED_ASSEMBLY_FILE="$POLISHED_ASSEMBLY_DIR/Bolinus_k63_polished.fasta"

# BUSCO lineage
BUSCO_LINEAGE="/data/proj2/home/students/m.borgmann/Master_thesis/data/raw/Quality_control/Assembly_stats/BUSCO_analysis/lineage_datasets/mollusca_odb10"

# Settings
THREADS=60

# Create directories if they don't exist
mkdir -p "$POLISHED_ASSEMBLY_DIR" "$MAPPING_DIR" "$QUAST_OUTPUT" "$BUSCO_OUTPUT"

# Function to check if a step needs to be run
check_step() {
    if [ -f "$1" ]; then
        echo "Found $1, skipping previous steps"
        return 0
    fi
    return 1
}

echo "Checking for polished assembly at $POLISHED_ASSEMBLY_FILE"
if [ -f "$POLISHED_ASSEMBLY_FILE" ]; then
    echo "Found polished assembly."
else
    echo "Polished assembly not found."

    # Check for sorted BAM file
    if ! check_step "$MAPPING_DIR/sorted_reads.bam"; then
        echo "Sorted BAM file not found."

        # Check for existing BAM file
        if check_step "$MAPPING_DIR/aligned_reads.bam"; then
            echo "Found BAM file at $MAPPING_DIR/aligned_reads.bam. Proceeding to sorting step."
            echo "Sorting BAM file..."
            samtools sort -@ "$THREADS" "$MAPPING_DIR/aligned_reads.bam" -o "$MAPPING_DIR/sorted_reads.bam" \
                > "$MAPPING_DIR/samtools_sort.log" 2>&1

            echo "Indexing sorted BAM file..."
            samtools index -@ "$THREADS" "$MAPPING_DIR/sorted_reads.bam" \
                > "$MAPPING_DIR/samtools_index.log" 2>&1

        elif check_step "$MAPPING_DIR/aligned_reads.sam"; then
            echo "Found SAM file at $MAPPING_DIR/aligned_reads.sam."
            echo "Converting SAM to BAM..."
            samtools view -@ "$THREADS" -Sb "$MAPPING_DIR/aligned_reads.sam" -o "$MAPPING_DIR/aligned_reads.bam" \
                > "$MAPPING_DIR/samtools_view.log" 2>&1

            echo "Sorting BAM file..."
            samtools sort -@ "$THREADS" "$MAPPING_DIR/aligned_reads.bam" -o "$MAPPING_DIR/sorted_reads.bam" \
                > "$MAPPING_DIR/samtools_sort.log" 2>&1

            echo "Indexing sorted BAM file..."
            samtools index -@ "$THREADS" "$MAPPING_DIR/sorted_reads.bam" \
                > "$MAPPING_DIR/samtools_index.log" 2>&1
        else
            echo "SAM file not found. Starting mapping reads..."
            echo "Indexing assembly with bwa-mem2..."
            bwa-mem2 index "$ASSEMBLY_FILE"

            echo "Aligning reads with bwa-mem2..."
            bwa-mem2 mem -t "$THREADS" "$ASSEMBLY_FILE" "$READ1" "$READ2" > "$MAPPING_DIR/aligned_reads.sam"
            echo "Alignment completed."

            echo "Converting SAM to BAM..."
            samtools view -@ "$THREADS" -Sb "$MAPPING_DIR/aligned_reads.sam" -o "$MAPPING_DIR/aligned_reads.bam" \
                > "$MAPPING_DIR/samtools_view.log" 2>&1

            echo "Sorting BAM file..."
            samtools sort -@ "$THREADS" "$MAPPING_DIR/aligned_reads.bam" -o "$MAPPING_DIR/sorted_reads.bam" \
                > "$MAPPING_DIR/samtools_sort.log" 2>&1

            echo "Indexing sorted BAM file..."
            samtools index -@ "$THREADS" "$MAPPING_DIR/sorted_reads.bam" \
                > "$MAPPING_DIR/samtools_index.log" 2>&1
        fi
    else
        echo "Found sorted BAM file at $MAPPING_DIR/sorted_reads.bam"
    fi

    echo "Running Pilon..."
    java -Xmx256g -jar /data/proj2/home/students/m.borgmann/software/pilon/pilon-1.24.jar \
        --genome "$ASSEMBLY_FILE" \
        --frags "$MAPPING_DIR/sorted_reads.bam" \
        --output "$POLISHED_ASSEMBLY_DIR/Bolinus_k63_polished" \
        > "$OUTPUT_DIR/pilon.log" 2>&1

    echo "Pilon polishing completed."
fi

# Run QC only if polished assembly exists
if [ -f "$POLISHED_ASSEMBLY_FILE" ]; then
    echo "Running QUAST..."
    quast.py "$POLISHED_ASSEMBLY_FILE" \
        --threads "$THREADS" \
        --labels "Bolinus_SOAP_k63_polished" \
        --eukaryote \
        --min-contig 1000 \
        --large \
        -o "$QUAST_OUTPUT" \
        > "$QUAST_OUTPUT/quast.log" 2>&1

    echo "QUAST analysis completed."

    echo "Running BUSCO..."
    busco -i "$POLISHED_ASSEMBLY_FILE" \
        -l "$BUSCO_LINEAGE" \
        -o "Bolinus_SOAP_k63_polished" \
        -m genome \
        --out_path "$BUSCO_OUTPUT" \
        --cpu "$THREADS" \
        > "$BUSCO_OUTPUT/busco.log" 2>&1

    echo "BUSCO analysis completed."
else
    echo "Polished assembly not found. Skipping QC steps."
fi

echo "Pipeline completed on $(date)."
