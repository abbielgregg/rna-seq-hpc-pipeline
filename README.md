# Reproducibility Assessment of RNA-seq Quantification Pipelines

An end-to-end bioinformatics and statistical framework evaluating the impact of quantification strategies on differential gene expression profiling. 

This repository evaluates how the choice of software pipelines influences biological conclusions when analyzing the same raw genetic data. It utilizes public data from a genetically engineered mouse model focusing on prostate cancer initiation driven by the Trp53 R270H mutation.

## Project Importance

When scientists study diseases like cancer using RNA sequencing, they rely on computational pipelines to convert raw biological data into a list of genes turned up or down by a disease. 

However, different algorithms use vastly different mathematical and biological assumptions. If a gene is flagged as highly critical by one software package but completely ignored by another, downstream medical conclusions or therapeutic targets could be compromised. This project serves as a transparency benchmark to separate genuine biological signals from software-specific artifacts.

## What Was Done

The analysis was performed on a dataset consisting of three genotype groups: Wild-Type, Heterozygous, and Homozygous. The core benchmark focused on the Homozygous vs Wild-Type comparison across two major computational steps:

### 1. Gene Quantification
Compared two fundamentally different mapping strategies:
* Alignment-based: HISAT2 and StringTie, which maps raw data precisely onto a physical reference genome.
* Pseudoalignment-based: Kallisto and tximport, an alignment-free approach that breaks data into tiny fragments to drastically speed up calculations.

### 2. Differential Expression Frameworks
Tested three industry-standard tools on the exact same quantified data to observe differences in statistical mathematics:
* DESeq2: Negative Binomial modeling with shrinkage estimation.
* edgeR: Negative Binomial modeling with empirical Bayes estimation.
* limma-voom: Linear modeling designed initially for microarrays, adapted with weight transformations.

### 3. Biological and Isoform Validation
* Pathway Analysis: Used clusterProfiler to check if varying gene lists pointed to the same functional biological processes.
* Transcript-Level Resolution: Explored individual gene alternate splicing via DESeq2 to capture hidden signals masked at the general gene level.

## Summary of Results

### Pipeline Quantification Agreement
* High Global Correlation: Changing the underlying mapping strategy still preserved a strong global correlation of r = 0.767 across all 17,159 analyzed genes.
* Overlapping Targets: The two quantification approaches shared a 58.2 percent overlap (428 genes) of core significant altered genes.

### Statistical Software Divergence
The choice of statistical tool drastically scaled the volume of discoveries under identical cutoff criteria:
* DESeq2 identified 734 genes (525 up, 209 down).
* edgeR identified 632 genes (421 up, 211 down).
* limma-voom was highly conservative, identifying only 8 genes (7 up, 1 down). 
* Concordance: All 8 genes flagged by limma-voom were successfully captured by both DESeq2 and edgeR, proving they represent the absolute highest-confidence biological signals.

### Biological Interpretations Shift
* Microenvironment Alterations: DESeq2 and edgeR primarily highlighted structural changes in the tumor microenvironment, dominating their profiles with muscle system and tissue processes.
* Direct Mutation Signals: Because limma-voom only captured a tiny pool of highly significant genes, it bypassed the structural muscle background and directly pinpointed p53-class mediation and cellular apoptosis pathways.
* Isoform Insights: High-resolution transcript exploration successfully unmasked isoform switching behaviors in cancer-associated genes like Hnrnpa1 that were completely hidden during traditional gene-level evaluations.

## Repository Structure

```text
rna-seq-hpc-pipeline/
├── data/
│   └── Metadata.txt
├── scripts/
│   ├── 01_pipeline_quantification.sh
│   ├── 02_kallisto_quant.sh
│   ├── 03_hisat2_deseq2.R
│   ├── 04_hisat2_edger.R
│   ├── 05_hisat2_limma_voom.R
│   ├── 06_kallisto_deseq2.R
│   ├── 07_kallisto_sleuth.R
│   ├── 08_pca_generation.R                 
│   ├── 09_pipeline_benchmarking_plots.R    
│   ├── 10_de_method_benchmarking.R         
│   ├── 11_go_enrichment_dotplots.R         
│   └── 12_isoform_expression_plots.R       
└── README.md
```

## Infrastructure and Dependencies
* High-Performance Computing: Deployed and run via an Imperial College Cluster HPC environment.
* Upstream Linux Shell Tools: hisat2/2.2.1, samtools/1.6, kallisto/0.51.1, python/2.7.11.
* Downstream R Statistical Packages: R/4.5.2, DESeq2, edgeR, limma, tximport, sleuth, clusterProfiler, UpSetR, ggplot2.


