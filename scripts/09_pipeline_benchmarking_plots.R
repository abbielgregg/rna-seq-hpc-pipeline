# Script for pipeline correlation and DEG concordance analysis

# Script 09: Cross-Pipeline Quantification Benchmarking & Intersections
# Metrics: Pearson Correlation (Concordance Mapping) & Manual Geometric Venns
# Comparison: HISAT2-StringTie vs Kallisto-tximport (HOM vs WT)

library(ggplot2)
library(dplyr)

# Load Core Structural Datasets
hisat2   <- read.csv("~/Desktop/Degree/NEWdata/DATA/DESeq2/DESeq2_HOM_vs_WT_ALL.csv")
kallisto <- read.csv("~/Desktop/Degree/NEWdata/DATA/kallisto_tximport/kallisto_HOM_vs_WT_all_genes.csv")

# Clean, Harmonise and Strip Suffix Variants
hisat2_clean <- hisat2 %>%
  rename(gene = X) %>%
  select(gene, baseMean_hisat2 = baseMean, log2FC_hisat2 = log2FoldChange, padj_hisat2 = padj)

kallisto_clean <- kallisto %>%
  mutate(gene = sub("\\..*", "", gene)) %>%
  select(gene, baseMean_kallisto = baseMean, log2FC_kallisto = log2FoldChange, padj_kallisto = padj)

# Combine on shared gene indicators and filter low-expressed variance shifts
merged <- inner_join(hisat2_clean, kallisto_clean, by = "gene") %>%
  filter(baseMean_hisat2 > 10) %>%
  filter(!is.na(log2FC_hisat2), !is.na(log2FC_kallisto), !is.na(padj_hisat2), !is.na(padj_kallisto))

# Classify Reproducibility Boundaries
merged <- merged %>%
  mutate(
    sig_hisat2   = padj_hisat2 < 0.05 & abs(log2FC_hisat2) > 1,
    sig_kallisto = padj_kallisto < 0.05 & abs(log2FC_kallisto) > 1,
    category     = case_when(
      sig_hisat2 & sig_kallisto   ~ "Shared DEG",
      sig_hisat2 & !sig_kallisto  ~ "HISAT2-StringTie only",
      !sig_hisat2 & sig_kallisto  ~ "Kallisto-tximport only",
      TRUE                        ~ "Not significant"
    ),
    category = factor(category, levels = c("Not significant", "Shared DEG", "HISAT2-StringTie only", "Kallisto-tximport only"))
  )

r_val <- round(cor(merged$log2FC_hisat2, merged$log2FC_kallisto, use = "complete.obs"), 3)
n_val <- nrow(merged)

# VISUALISATION A: Concordance Scatter Mapping

scatter_plot <- ggplot(merged, aes(x = log2FC_hisat2, y = log2FC_kallisto, colour = category)) +
  geom_point(data = filter(merged, category == "Not significant"), alpha = 0.2, size = 0.6) +
  geom_point(data = filter(merged, category != "Not significant"), alpha = 0.8, size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.6) +
  scale_colour_manual(
    values = c(
      "Not significant"        = "grey60",
      "Shared DEG"             = "#0072B2", # Color-blind safe blue
      "HISAT2-StringTie only"  = "#E69F00", # Color-blind safe orange
      "Kallisto-tximport only" = "#009E73"  # Color-blind safe green
    ),
    drop = FALSE
  ) +
  coord_cartesian(xlim = c(-6, 6), ylim = c(-6, 6)) +
  annotate("text", x = -5.8, y = 5.8, label = paste0("r = ", r_val, "\nn = ", format(n_val, big.mark = ",")), 
           hjust = 0, vjust = 1, size = 3.8, fontface = "italic", colour = "black") +
  labs(x = "HISAT2-StringTie log2FC", y = "Kallisto-tximport log2FC", colour = NULL) +
  theme_classic(base_size = 12) +
  theme(
    legend.position   = "bottom",
    legend.text       = element_text(size = 10),
    axis.title        = element_text(size = 11, face = "bold"),
    axis.text         = element_text(size = 10),
    panel.background  = element_rect(fill = "white", colour = NA),
    plot.background   = element_rect(fill = "white", colour = NA),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank()
  ) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1), nrow = 2))

ggsave("~/Desktop/Degree/NEWdata/scatter_pipeline_comparison_REPORT.png", plot = scatter_plot, width = 5.5, height = 5.5, dpi = 300)

# VISUALISATION B: Manual Geometric Venn Diagram Integration

counts_table  <- table(merged$category)
hisat2_only   <- as.numeric(counts_table["HISAT2-StringTie only"])
kallisto_only <- as.numeric(counts_table["Kallisto-tximport only"])
shared        <- as.numeric(counts_table["Shared DEG"])
total_deg     <- hisat2_only + shared + kallisto_only

circle_radius <- 1.2
center_offset <- 0.75

venn_plot <- ggplot() +
  annotate("polygon", x = center_offset * -1 + circle_radius * cos(seq(0, 2*pi, length.out = 300)),
           y = circle_radius * sin(seq(0, 2*pi, length.out = 300)), fill = "#E69F00", alpha = 0.25) +
  annotate("polygon", x = center_offset * 1 + circle_radius * cos(seq(0, 2*pi, length.out = 300)),
           y = circle_radius * sin(seq(0, 2*pi, length.out = 300)), fill = "#009E73", alpha = 0.25) +
  annotate("path", x = center_offset * -1 + circle_radius * cos(seq(0, 2*pi, length.out = 300)),
           y = circle_radius * sin(seq(0, 2*pi, length.out = 300)), colour = "#E69F00", linewidth = 1.2, alpha = 0.9) +
  annotate("path", x = center_offset * 1 + circle_radius * cos(seq(0, 2*pi, length.out = 300)),
           y = circle_radius * sin(seq(0, 2*pi, length.out = 300)), colour = "#009E73", linewidth = 1.2, alpha = 0.9) +
  annotate("text", x = -1.1, y = 0.15, label = hisat2_only, size = 7, fontface = "bold", colour = "black", hjust = 0.5) +
  annotate("text", x = -1.1, y = -0.15, label = paste0("(", round(hisat2_only/total_deg*100, 1), "%)"), size = 3.8, colour = "black", hjust = 0.5) +
  annotate("text", x = 0, y = 0.15, label = shared, size = 7, fontface = "bold", colour = "black", hjust = 0.5) +
  annotate("text", x = 0, y = -0.15, label = paste0("(", round(shared/total_deg*100, 1), "%)"), size = 3.8, colour = "black", hjust = 0.5) +
  annotate("text", x = 1.1, y = 0.15, label = kallisto_only, size = 7, fontface = "bold", colour = "black", hjust = 0.5) +
  annotate("text", x = 1.1, y = -0.2, label = paste0("(", round(kallisto_only/total_deg*100, 1), "%)"), size = 3.8, colour = "black", hjust = 0.5) +
  annotate("text", x = -0.75, y = 1.5, label = "HISAT2-StringTie", size = 4, fontface = "bold", colour = "black", hjust = 0.5) +
  annotate("text", x = 0.75, y = 1.5, label = "Kallisto-tximport", size = 4, fontface = "bold", colour = "black", hjust = 0.5) +
  theme_void() + coord_fixed()

ggsave("~/Desktop/Degree/NEWdata/Venn_pipeline_comparison.png", plot = venn_plot, width = 5, height = 4, dpi = 300)
cat("Benchmarking visualization elements compiled completely!\n")

