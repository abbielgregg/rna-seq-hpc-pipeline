# Script for statistical method comparisons and sensitivity analysis

library(UpSetR)
library(ggplot2)
library(tidyr)
library(dplyr)

# FIGURE 4A: UpSet Plot
# Create a named list of your significant DEG vectors (FDR < 0.05, |log2FC| > 1)
deg_lists <- list(
  DESeq2 = deseq2_deg_vector,    # Replace with your actual gene list vector
  edgeR = edger_deg_vector,      # Replace with your actual gene list vector
  limma = limma_deg_vector       # Replace with your actual gene list vector
)

png("plots/Figure4A_UpSet_DE_Comparison.png", width = 600, height = 400)
upset(fromList(deg_lists), 
      order.by = "freq", 
      main.bar.color = "black", 
      sets.bar.color = "black",
      matrix.color = "black",
      shade.color = "grey80",
      text.scale = 1.3)
dev.off()


# FIGURE 4B: Directional Stacked Bar Chart
# Build a summary dataframe matching your counts (734, 632, 8)
direction_data <- data.frame(
  Method = c("DESeq2", "DESeq2", "edgeR", "edgeR", "limma", "limma"),
  Direction = c("Upregulated", "Downregulated", "Upregulated", "Downregulated", "Upregulated", "Downregulated"),
  Count = c(525, 209, 421, 211, 7, 1) # Populated from your report page 7
)

ggplot(direction_data, aes(x = Method, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", width = 0.5, color = "black") +
  scale_fill_manual(values = c("Upregulated" = "#3182bd", "Downregulated" = "#e6550d")) + # Example colors
  theme_classic() +
  labs(y = "Number of DEGs", x = "DE Method") +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14))

ggsave("plots/Figure4B_DEG_Direction.png", width = 5, height = 5)


# FIGURE 4C: FDR Sensitivity Analysis Line Plot
# Dataframe layout for tracking pipeline thresholds
fdr_sensitivity <- data.frame(
  FDR_Threshold = rep(c(0.001, 0.005, 0.01, 0.05, 0.1), 3),
  Method = rep(c("DESeq2", "edgeR", "limma-voom"), each = 5),
  Percent_DE = c(
    # DESeq2 % values across 0.001 to 0.1
    # edgeR % values across 0.001 to 0.1
    # limma-voom % values across 0.001 to 0.1
  )
)

ggplot(fdr_sensitivity, aes(x = FDR_Threshold, y = Percent_DE, color = Method, group = Method)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_x_log10(breaks = c(0.001, 0.005, 0.01, 0.05, 0.1)) + 
  theme_classic() +
  labs(x = "FDR threshold", y = "Genes classified as DE (%)", color = "Method")

