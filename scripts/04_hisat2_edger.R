# Script 04: Differential Expression Analysis via edgeR (HISAT2 Count Matrix)
# Statistical Framework: Empirical Bayes Dispersion Partitioning models

library(edgeR)

counts <- read.csv("gene_count_matrix.csv", row.names=1)
meta   <- read.table("data/Metadata.txt", header=TRUE, sep="\t", row.names=1)

counts <- counts[, rownames(meta)]
group  <- factor(meta$condition, levels=c("WT","HET","HOM"))

# Build core edgeR collection environment container
y <- DGEList(counts=counts, group=group)

# Apply internal expression scaling metrics
keep <- filterByExpr(y)
y    <- y[keep, , keep.lib.sizes=FALSE]
y    <- calcNormFactors(y)

# Construct design configuration matrix
design <- model.matrix(~group)
y      <- estimateDisp(y, design)
fit    <- glmFit(y, design)

# Evaluate target coefficient conditions
lrt_HOM_vs_WT  <- glmLRT(fit, coef="groupHOM")
lrt_HET_vs_WT  <- glmLRT(fit, coef="groupHET")
lrt_HOM_vs_HET <- glmLRT(fit, contrast=c(0, -1, 1))

# Unpack baseline tables
res_HOM_vs_WT  <- topTags(lrt_HOM_vs_WT, n=Inf)$table
res_HET_vs_WT  <- topTags(lrt_HET_vs_WT, n=Inf)$table
res_HOM_vs_HET <- topTags(lrt_HOM_vs_HET, n=Inf)$table

write.csv(res_HOM_vs_WT,  "edgeR_HOM_vs_WT.csv")
write.csv(res_HET_vs_WT,  "edgeR_HET_vs_WT.csv")
write.csv(res_HOM_vs_HET, "edgeR_HOM_vs_HET.csv")

# Filter significantly altered gene matrices (FDR < 0.05 and |logFC| > 1)
filter_edgeR <- function(df) {
  df[df$FDR < 0.05 & abs(df$logFC) > 1, ]
}

write.csv(filter_edgeR(res_HOM_vs_WT),  "edgeR_HOM_vs_WT_significant.csv")
write.csv(filter_edgeR(res_HET_vs_WT),  "edgeR_HET_vs_WT_significant.csv")
write.csv(filter_edgeR(res_HOM_vs_HET), "edgeR_HOM_vs_HET_significant.csv")

