# Pipeline execution script for raw read alignment using HISAT2

#!/bin/bash
# Primary scripts folder containing core quantification and downstream statistical analysis pipelines
# Script 01: Preprocessing, Alignment and Alignment-Based Expression Quantification
# Workflow: HISAT2 (Genome Mapping) -> Samtools (Sorting/Indexing) -> StringTie -> Python Matrix

# Server Environment Setup (Executed via SSH on Imperial HPC)
# Connection String: ssh ag2022@compnode.bc.ic.ac.uk
module load hisat2/2.2.1
module load samtools/1.6
module load python/2.7.11

# Directory Path Nicknames
BASE_DIR="/users/lux/lqp17/rna-seq-project"
DATA_DIR="/project/data/huntley/lh-project/ma-rg/reads"
ALIGN_DIR="${BASE_DIR}/alignments"
STRINGTIE_DIR="${BASE_DIR}/stringtie"

mkdir -p ${ALIGN_DIR}
mkdir -p ${STRINGTIE_DIR}

# Reference Genome Paths (GRCm39 Primary Assembly Release)
GENOME_INDEX="/path/to/hisat2-index/mouse-grcm39"
GTF_ANNOTATION="/path/to/annotation/GRCm39.113.gtf"

# Experimental Replicates Target List
SAMPLES=("SRR8925047" "SRR8925048" "SRR8925049" "SRR8925050" "SRR8925051" "SRR8925052")

# Core Automation Alignment Loop
# Managed via background screening to survive dropouts: screen -S rna_seq_run
for SAMPLE in "${SAMPLES[@]}"
do
    # Align Paired-End Fastq Reads directly into sorted BAM binaries to minimize storage footprints
    hisat2 -x ${GENOME_INDEX} \
           -1 ${DATA_DIR}/${SAMPLE}_1.fastq.gz \
           -2 ${DATA_DIR}/${SAMPLE}_2.fastq.gz \
           | samtools view -bS - \
           | samtools sort -o ${ALIGN_DIR}/${SAMPLE}_sorted.bam -
           
    # Index the sorted file for quick random spatial reads access
    samtools index ${ALIGN_DIR}/${SAMPLE}_sorted.bam
    
    # Assemble structural transcripts and capture sample calculations
    mkdir -p ${STRINGTIE_DIR}/${SAMPLE}
    stringtie ${ALIGN_DIR}/${SAMPLE}_sorted.bam \
              -e -B \
              -G ${GTF_ANNOTATION} \
              -o ${STRINGTIE_DIR}/${SAMPLE}/${SAMPLE}_merged.gtf
done

# Consolidate Results into Structural Master Matrices
cat <<EOF > ${BASE_DIR}/samples.txt
SRR8925047 ${STRINGTIE_DIR}/SRR8925047/SRR8925047_merged.gtf
SRR8925048 ${STRINGTIE_DIR}/SRR8925048/SRR8925048_merged.gtf
SRR8925049 ${STRINGTIE_DIR}/SRR8925049/SRR8925049_merged.gtf
SRR8925050 ${STRINGTIE_DIR}/SRR8925050/SRR8925050_merged.gtf
SRR8925051 ${STRINGTIE_DIR}/SRR8925051/SRR8925051_merged.gtf
SRR8925052 ${STRINGTIE_DIR}/SRR8925052/SRR8925052_merged.gtf
EOF

python prepDE.py -i ${BASE_DIR}/samples.txt \
                 -g gene_count_matrix.csv \
                 -t transcript_count_matrix.csv

