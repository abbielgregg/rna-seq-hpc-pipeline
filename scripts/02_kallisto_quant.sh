#!/bin/bash
# Script 02: Alignment-Free Transcript Pseudoalignment & Inference Estimation
# Workflow: Kallisto Quant (100 Bootstrap Iterations for Uncertainty Resolution)

module load kallisto/derek

BASE_DIR="/users/y3/ag2622/kallisto_tximport"
READS_DIR="/project/data/huntley/fy_projects/rna-seq/reads"
INDEX_PATH="/project/data/huntley/fy_projects/rna-seq/kallisto_index/index.idx"
OUTPUT_DIR="${BASE_DIR}/quant"

mkdir -p ${OUTPUT_DIR}

SAMPLES=("SRR8985047" "SRR8985048" "SRR8985049" "SRR8985050" "SRR8985051" "SRR8985052")

for s in "${SAMPLES[@]}"
do
    mkdir -p ${OUTPUT_DIR}/${s}
    
    # Run pseudoalignment tracking with 100 inferential boot-replicates for sleuth functionality
    kallisto quant \
        -i ${INDEX_PATH} \
        -o ${OUTPUT_DIR}/${s} \
        -b 100 \
        ${READS_DIR}/${s}_1.fastq.gz \
        ${READS_DIR}/${s}_2.fastq.gz
done

