library(ggplot2)

# Build manual evaluation frame from your transcript results summary table
isoform_df <- data.frame(
  Transcript = c("ENSMUST00000253033", "ENSMUST00000253035", "ENSMUST00000253039", "ENSMUST00000253034"),
  log2FC = c(-2.0, -2.1, -2.4, -2.8), # Add exact data points from your report source matrix
  SE = c(0.3, 0.25, 0.35, 0.4),       # Add standard errors calculated from results
  Direction = "Down"
)

ggplot(isoform_df, aes(x = log2FC, y = reorder(Transcript, log2FC), fill = Direction)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  geom_errorbar(aes(xmin = log2FC - SE, xmax = log2FC + SE), width = 0.2) +
  scale_fill_manual(values = c("Down" = "#df9f28")) + # Matching report palette
  theme_classic() +
  labs(x = "log2 Fold Change (HOM vs WT)", y = "Transcript Isoform", title = "Isoform-level expression of ENSMUSG00000122428") +
  theme(legend.position = "right")

