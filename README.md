# Reproducibility Assessment of RNA-seq Quantification Pipelines

An end-to-end bioinformatics and statistical framework evaluating the impact of quantification strategies on differential gene expression (DGE) profiling. This project utilizes paired-end RNA-seq data from mouse models (*Mus musculus* GRCm39) carrying the **Trp53 R270H mutation** across Wild-Type (WT), Heterozygous (HET), and Homozygous (HOM) cohorts to model the boundaries of pipeline sensitivity.

##  Core Objectives
*   **Pipeline Benchmark:** Direct cross-method evaluation comparing standard alignment-based pipelines (**HISAT2 + StringTie**) against high-throughput alignment-free pseudoalignment tools (**Kallisto**).
*   **Downstream Sensitivity:** Assessment of statistical variations across standard differential expression frameworks (**DESeq2**, **edgeR**, and **limma-voom**).
*   **Resolution Modeling:** Isolating transcript-isoform expression mechanics using uncertainty-aware bootstrap variance metrics (**Sleuth**).

---

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
│   ├── 08_pca_generation.R              # <-- NEW: Direct VST-PCA Mapping Plotter
│   └── 09_pipeline_benchmarking_plots.R # <-- NEW: Scatter Concordance & Venn Plots
└── README.md
```

## ️ Infrastructure & Dependencies
*   **High-Performance Computing:** Deployed via Imperial College Cluster HPC environment.
*   **Upstream Core Linux Shell Tools:** `hisat2/2.2.1`, `samtools/1.6`, `kallisto/derek`, `python/2.7.11`
*   **Downstream R Statistical Packages:** `DESeq2`, `edgeR`, `limma`, `tximport`, `sleuth`

---

