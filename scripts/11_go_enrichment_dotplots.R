# Script for functional pathway enrichment analysis using clusterProfiler

library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)

# Run this block uniquely for your DESeq2, edgeR, limma (gene), and DESeq2 (transcript) vectors
go_enrich <- enrichGO(
  gene          = deseq2_deg_entrez_ids, # Needs Entrez IDs as mapped by org.Mm.eg.db
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)

# Plot standard dotplot matching layout
dotplot(go_enrich, showCategory = 15) +
  scale_color_gradient(low = "red", high = "blue", name = "p.adjust") + 
  theme_classic() +
  labs(title = "GO BP enrichment - DESeq2 (HOM vs WT)")

