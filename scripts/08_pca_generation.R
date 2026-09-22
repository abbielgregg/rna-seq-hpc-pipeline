# Script for sample quality control and PCA visualization

# Script 08: Comparative Principal Component Analysis (PCA) Generation
# Workflow: DESeq2 VST Normalisation -> ggplot2 Publication Visualisation
# Comparison: Alignment-Based (HISAT2) vs Alignment-Free (Kallisto) Layouts

library(DESeq2)
library(tximport)
library(ggplot2)

# Master Plot Aesthetics & Comparable Boundaries
color_palette <- c("WT" = "#44AA99", "HET" = "#56B4E9", "HOM" = "#D55E00")
shape_palette <- c("WT" = 16, "HET" = 17, "HOM" = 15)
fixed_x_window <- c(-40, 40)
fixed_y_window <- c(-20, 20)

shared_pca_theme <- theme_classic(base_size = 13) +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0.5),
    axis.title   = element_text(face = "bold"),
    axis.text    = element_text(colour = "black"),
    legend.title = element_text(face = "bold")
  )

# PANEL A: HISAT2-StringTie PCA

counts_hisat2      <- read.csv("~/Desktop/Degree/NEWdata/gene_count_matrix.csv", row.names = 1)
sample_info_hisat2 <- data.frame(
  row.names = c("SRR8985047", "SRR8985048", "SRR8985049", "SRR8985050", "SRR8985051", "SRR8985052"),
  condition = factor(c("WT", "WT", "HET", "HET", "HOM", "HOM"), levels = c("WT", "HET", "HOM"))
)

dds_hisat2 <- DESeqDataSetFromMatrix(countData = round(counts_hisat2), colData = sample_info_hisat2, design = ~ condition)
dds_hisat2 <- dds_hisat2[rowSums(counts(dds_hisat2)) >= 10, ] # Filter background noise
vst_hisat2 <- vst(dds_hisat2, blind = TRUE)

pca_data_hisat2    <- plotPCA(vst_hisat2, intgroup = "condition", returnData = TRUE)
percentVar_hisat2  <- round(100 * attr(pca_data_hisat2, "percentVar"))

plot_hisat2_pca <- ggplot(pca_data_hisat2, aes(x = PC1, y = PC2, colour = condition, shape = condition)) +
  geom_point(size = 4.5, alpha = 0.9) +
  scale_colour_manual(values = color_palette) +
  scale_shape_manual(values = shape_palette) +
  coord_cartesian(xlim = fixed_x_window, ylim = fixed_y_window) +
  labs(
    title   = "HISAT2-StringTie",
    x       = paste0("PC1: ", percentVar_hisat2[1], "% variance"),
    y       = paste0("PC2: ", percentVar_hisat2[2], "% variance"),
    colour  = "Condition",
    shape   = "Condition"
  ) +
  shared_pca_theme

ggsave("~/Desktop/Degree/NEWdata/PCA_hisat2.png", plot = plot_hisat2_pca, width = 5, height = 5, dpi = 300)

# PANEL B: Kallisto-tximport PCA
tx2gene <- read.table("~/Desktop/Degree/NEWdata/t2g.txt", header = FALSE, col.names = c("gene_id", "transcript_id"))
tx2gene <- tx2gene[, c("transcript_id", "gene_id")] # Force proper tximport index tracking

sample_names <- c("SRR8985047", "SRR8985048", "SRR8985049", "SRR8985050", "SRR8985051", "SRR8985052")
files        <- file.path("~/Desktop/Degree/NEWdata/quant", sample_names, "abundance.h5")
names(files) <- sample_names

txi <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion = TRUE)

sample_info_kallisto <- data.frame(
  row.names = c("SRR8985047", "SRR8985048", "SRR8985049", "SRR8985050", "SRR8985051", "SRR8985052"),
  condition = factor(c("WT", "WT", "HET", "HET", "HOM", "HOM"), levels = c("WT", "HET", "HOM"))
)

dds_kallisto <- DESeqDataSetFromTximport(txi = txi, colData = sample_info_kallisto, design = ~ condition)
dds_kallisto <- dds_kallisto[rowSums(counts(dds_kallisto)) >= 10, ]
vst_kallisto <- vst(dds_kallisto, blind = TRUE)

pca_data_kallisto   <- plotPCA(vst_kallisto, intgroup = "condition", returnData = TRUE)
percentVar_kallisto <- round(100 * attr(pca_data_kallisto, "percentVar"))

plot_kallisto_pca <- ggplot(pca_data_kallisto, aes(x = PC1, y = PC2, colour = condition, shape = condition)) +
  geom_point(size = 4.5, alpha = 0.9) +
  scale_colour_manual(values = color_palette) +
  scale_shape_manual(values = shape_palette) +
  coord_cartesian(xlim = fixed_x_window, ylim = fixed_y_window) +
  labs(
    title   = "Kallisto-tximport",
    x       = paste0("PC1: ", percentVar_kallisto[1], "% variance"),
    y       = paste0("PC2: ", percentVar_kallisto[2], "% variance"),
    colour  = "Condition",
    shape   = "Condition"
  ) +
  shared_pca_theme

ggsave("~/Desktop/Degree/NEWdata/PCA_kallisto.png", plot = plot_kallisto_pca, width = 5, height = 5, dpi = 300)
cat("PCA generation completed! Files saved to output paths.\n")

