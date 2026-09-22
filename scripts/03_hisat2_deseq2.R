# Script 03: Differential Expression Analysis via DESeq2 (HISAT2 Count Matrix)
# Statistical Framework: Negative Binomial Generalized Linear Models

library(DESeq2)

# Load cross-matched parameters
counts <- read.csv("gene_count_matrix.csv", row.names=1)
meta   <- read.table("data/Metadata.txt", header=TRUE, sep="\t", row.names=1)

# Sort coordinates alignment check
counts <- counts[, rownames(meta)]
meta$condition <- factor(meta$condition, levels=c("WT","HET","HOM"))

# Build design construct
dds <- DESeqDataSetFromMatrix(countData=counts, colData=meta, design=~condition)

# Filter low count noise to reduce multi-testing FDR penalties
dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)

# Process cross-comparisons results matrices
res_HOM_vs_WT  <- results(dds, contrast=c("condition","HOM","WT"))
res_HET_vs_WT  <- results(dds, contrast=c("condition","HET","WT"))
res_HOM_vs_HET <- results(dds, contrast=c("condition","HOM","HET"))

# Save master unfiltered datasets
write.csv(as.data.frame(res_HOM_vs_WT),  "DESeq2_HOM_vs_WT_ALL.csv")
write.csv(as.data.frame(res_HET_vs_WT),  "DESeq2_HET_vs_WT_ALL.csv")
write.csv(as.data.frame(res_HOM_vs_HET), "DESeq2_HOM_vs_HET_ALL.csv")

# Filter statistically robust targets (|log2FC| > 1 and adjusted p-value < 0.05)
filter_sig <- function(res) {
  res[!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) > 1, ]
}

write.csv(as.data.frame(filter_sig(res_HOM_vs_WT)),  "DESeq2_HOM_vs_WT_sig.csv")
write.csv(as.data.frame(filter_sig(res_HET_vs_WT)),  "DESeq2_HET_vs_WT_sig.csv")
write.csv(as.data.frame(filter_sig(res_HOM_vs_HET)), "DESeq2_HOM_vs_HET_sig.csv")

# Persist active dds workspace instance for figure integration steps later
saveRDS(dds, "dds_hisat2_master.rds")

