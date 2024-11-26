#!/bin/bash

# Paths
TRIMMOMATIC_JAR=/data/proj2/home/students/m.borgmann/software/Trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar
ADAPTERS=/data/proj2/home/students/m.borgmann/software/Trimmomatic/Trimmomatic-0.39/adapters/TruSeq3-PE.fa
FASTP=/data/proj2/home/students/m.borgmann/miniforge3/envs/master_thesis/bin/fastp
INPUT_DIR=/data/proj2/home/students/m.borgmann/Master_thesis/data/raw/Reference_data/WGS_projects
OUTPUT_DIR=/data/proj2/home/students/m.borgmann/Master_thesis/data/processed/QCed_IlluminaFiles

# Trimming parameters
LEADING=3
TRAILING=3
SLIDINGWINDOW=4:15
MINLEN=36

# Function to process a single sample
process_sample() {
    local PROJECT=$1
    local SAMPLE=$2
    local SPECIES_DIR=$3

    # Input and output file paths for Trimmomatic
    local INPUT1=${INPUT_DIR}/${PROJECT}/${SAMPLE}/${SAMPLE}_1.fastq
    local INPUT2=${INPUT_DIR}/${PROJECT}/${SAMPLE}/${SAMPLE}_2.fastq
    local TRIMMED1=${SPECIES_DIR}/${SAMPLE}_1.trimmed.fastq
    local TRIMMED2=${SPECIES_DIR}/${SAMPLE}_2.trimmed.fastq
    
    # Output files for fastp
    local FINAL1=${SPECIES_DIR}/${SAMPLE}_1.qcprocessed.fastq
    local FINAL2=${SPECIES_DIR}/${SAMPLE}_2.qcprocessed.fastq

    # Check if final FastQC reports exist
    if [[ -f ${SPECIES_DIR}/fastqc_postprocess_reports/${SAMPLE}_1.qcprocessed_fastqc.html && -f ${SPECIES_DIR}/fastqc_postprocess_reports/${SAMPLE}_2.qcprocessed_fastqc.html ]]; then
        echo "Post-processing FastQC reports already exist for sample: $SAMPLE. Skipping all processing."
        return 0
    fi

    # Check if trimmed files exist - if not, run Trimmomatic
    if [[ ! -f $TRIMMED1 || ! -f $TRIMMED2 ]]; then
        echo "Trimmed files not found. Processing sample with Trimmomatic: $SAMPLE from project: $PROJECT"
        
        java -jar $TRIMMOMATIC_JAR PE \
        -phred33 \
        -threads 30 \
        $INPUT1 $INPUT2 \
        $TRIMMED1 /dev/null \
        $TRIMMED2 /dev/null \
        ILLUMINACLIP:$ADAPTERS:2:30:10 \
        LEADING:$LEADING TRAILING:$TRAILING \
        SLIDINGWINDOW:$SLIDINGWINDOW \
        MINLEN:$MINLEN

        if [[ ! -f $TRIMMED1 || ! -f $TRIMMED2 ]]; then
            echo "Trimming failed for $SAMPLE. Please check input files and parameters."
            return 1
        fi
        
        echo "Running FastQC on trimmed files for sample: $SAMPLE"
        mkdir -p ${SPECIES_DIR}/fastqc_reports
        fastqc -t 30 -o ${SPECIES_DIR}/fastqc_reports $TRIMMED1 $TRIMMED2
    else
        echo "Trimmed files already exist for sample: $SAMPLE. Using these for polyG removal."
    fi

    # Run fastp on the trimmed files
    if [[ ! -f $FINAL1 || ! -f $FINAL2 ]]; then
        echo "Processing polyG removal for sample: $SAMPLE"
        mkdir -p ${SPECIES_DIR}/fastp_reports

        $FASTP \
        -i $TRIMMED1 \
        -I $TRIMMED2 \
        -o $FINAL1 \
        -O $FINAL2 \
        --trim_poly_g \
        --disable_adapter_trimming \
        --disable_quality_filtering \
        --disable_length_filtering \
        --thread 30 \
        --html ${SPECIES_DIR}/fastp_reports/${SAMPLE}.fastp_report.html \
        --json ${SPECIES_DIR}/fastp_reports/${SAMPLE}.fastp_report.json

        if [[ ! -f $FINAL1 || ! -f $FINAL2 ]]; then
            echo "PolyG removal failed for $SAMPLE. Please check trimmed files."
            return 1
        fi
    else
        echo "PolyG-trimmed files already exist for sample: $SAMPLE."
    fi

    # Run final FastQC on polyG-trimmed files
    echo "Running FastQC on polyG-trimmed files for sample: $SAMPLE"
    mkdir -p ${SPECIES_DIR}/fastqc_postprocess_reports
    fastqc -t 30 -o ${SPECIES_DIR}/fastqc_postprocess_reports $FINAL1 $FINAL2

    if [[ ! -f ${SPECIES_DIR}/fastqc_postprocess_reports/${SAMPLE}_1.qcprocessed_fastqc.html || ! -f ${SPECIES_DIR}/fastqc_postprocess_reports/${SAMPLE}_2.qcprocessed_fastqc.html ]]; then
        echo "Final FastQC failed for $SAMPLE. Please check qcprocessed files."
        return 1
    fi

    echo "Finished processing sample: $SAMPLE"
    return 0
}

# Create main output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Process Hexaplex trunculus samples
HEXAPLEX_DIR=${OUTPUT_DIR}/Hexaplex
mkdir -p $HEXAPLEX_DIR
process_sample "PRJNA1106542_Hexaplex_trunculus" "SRR28865916" "$HEXAPLEX_DIR"

# Process Bolinus brandaris samples
BOLINUS_DIR=${OUTPUT_DIR}/Bolinus
mkdir -p $BOLINUS_DIR
process_sample "PRJNA1106534_Bolinus_brandaris" "SRR28863561" "$BOLINUS_DIR"

echo "All processing completed. Check ${OUTPUT_DIR} for final QC processed files and reports."