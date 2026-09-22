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
│   └── Metadata.txt                 # Sample ID condition mapping configurations
├── scripts/
│   ├── 01_pipeline_quantification.sh # HISAT2 mapping & StringTie count aggregation
│   ├── 02_kallisto_quant.sh          # Kallisto pseudoalignment with 100 bootstrap layers
│   ├── 03_hisat2_deseq2.R            # DESeq2 computations on alignment metrics
│   ├── 04_hisat2_edger.R             # edgeR calculations on alignment matrices
│   ├── 05_hisat2_limma_voom.R        # Limma-voom evaluation across small cohorts
│   ├── 06_kallisto_deseq2.R          # tximport gene-level cohesion + DESeq2 modeling
│   └── 07_kallisto_sleuth.R          # Isoform-level variance tracking via Sleuth Wald tests
└── README.md                        # Project documentation and architecture logs
```

## ️ Infrastructure & Dependencies
*   **High-Performance Computing:** Deployed via Imperial College Cluster HPC environment.
*   **Upstream Core Linux Shell Tools:** `hisat2/2.2.1`, `samtools/1.6`, `kallisto/derek`, `python/2.7.11`
*   **Downstream R Statistical Packages:** `DESeq2`, `edgeR`, `limma`, `tximport`, `sleuth`

---

## Current Progress
1.  **Upstream Processing [COMPLETED]:** Achieved ~78-81% mapping across read cohorts. Generated master genomic count matrices via `prepDE.py` and abundance tracking variables via `.h5` containers.
2.  **Statistical Profiling [COMPLETED]:** Modeled 3-way physiological contrasts (HOM vs WT, HET vs WT, HOM vs HET) across all variations of negative-binomial, empirical Bayes, and linear weight distribution engines.
3.  **Visualization & Annotation Synthesis [UPCOMING]:** Cross-mapping Ensembl stable IDs to functional gene symbols and generating comprehensive pipeline overlap analyses (UpSet plots, Volcano metrics, PCA spacing).

