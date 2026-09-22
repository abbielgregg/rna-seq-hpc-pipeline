# Script 05: Differential Expression Analysis via Limma-Voom (HISAT2 Count Matrix)
# Statistical Framework: Linear Modeling with Variance Weights Allocation

library(limma)
library(edgeR)

counts <- read.csv("gene_count_matrix.csv", row.names=1)
meta   <- read.table("data/Metadata.txt", header=TRUE, sep="\t", row.names=1)

counts <- counts[, rownames(meta)]
group  <- factor(meta$condition, levels=c("WT","HET","HOM"))

y    <- DGEList(counts=counts, group=group)
keep <- filterByExpr(y)
y    <- y[keep, , keep.lib.sizes=FALSE]
y    <- calcNormFactors(y)

# Setup specialized 0-intercept design matrix formula
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# Transform log-counts per million measurements via voom weighting trends
v   <- voom(y, design, plot=FALSE)
fit <- lmFit(v, design)

# Map specific comparison linear vectors
contrasts <- makeContrasts(
  HOM_vs_WT  = HOM - WT,
  HET_vs_WT  = HET - WT,
  HOM_vs_HET = HOM - HET,
  levels     = design
)

fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)

# Pull full linear contrast summaries
res_HOM_vs_WT  <- topTable(fit2, coef="HOM_vs_WT", number=Inf)
res_HET_vs_WT  <- topTable(fit2, coef="HET_vs_WT", number=Inf)
res_HOM_vs_HET <- topTable(fit2, coef="HOM_vs_HET", number=Inf)

# Isolate passing significance boundaries (adj.P.Val < 0.05 and |logFC| > 1)
filter_limma <- function(df) {
  df[df$adj.P.Val < 0.05 & abs(df$logFC) > 1, ]
}

write.csv(filter_limma(res_HOM_vs_WT),  "limma_HOM_vs_WT_significant.csv")
write.csv(filter_limma(res_HET_vs_WT),  "limma_HET_vs_WT_significant.csv")
write.csv(filter_limma(res_HOM_vs_HET), "limma_HOM_vs_HET_significant.csv")

