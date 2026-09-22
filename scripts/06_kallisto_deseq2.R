# Script 06: Gene-Level Analysis via Kallisto + DESeq2
# Statistical Framework: tximport Transcript Cohesion Transformation

library(tximport)
library(readr)
library(DESeq2)
library(dplyr)
library(tibble)

base_dir    <- "/users/y3/ag2622/kallisto_tximport"
quant_dir   <- file.path(base_dir, "quant")
results_dir <- file.path(base_dir, "results")

samples           <- read.csv(file.path(base_dir, "samples.csv"), stringsAsFactors = FALSE)
samples$condition <- factor(samples$condition, levels = c("WT", "HET", "HOM"))

files        <- file.path(quant_dir, samples$sample, "abundance.h5")
names(files) <- samples$sample

# Ingest structural transcript index guides
tx2gene_full      <- read.delim("/project/data/huntley/fy_projects/rna-seq/kallisto_index/t2g.txt", header = FALSE, sep = "\t")
tx2gene           <- tx2gene_full[, 1:2]
colnames(tx2gene) <- c("TXNAME", "GENEID")
tx2gene           <- unique(tx2gene[!is.na(tx2gene$TXNAME) & !is.na(tx2gene$GENEID), ])

# Trim variant version tags to establish safe lookup matches
tx2gene$TXNAME <- sub("\\..*$", "", tx2gene$TXNAME)
tx2gene$GENEID <- sub("\\..*$", "", tx2gene$GENEID)
tx2gene        <- unique(tx2gene)

# Summarize abundances seamlessly to gene boundaries
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion = TRUE)

coldata <- samples %>% column_to_rownames("sample")
dds     <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ condition)
dds     <- DESeq(dds)

res_HET_vs_WT  <- results(dds, contrast = c("condition", "HET", "WT"))
res_HOM_vs_WT  <- results(dds, contrast = c("condition", "HOM", "WT"))
res_HOM_vs_HET <- results(dds, contrast = c("condition", "HOM", "HET"))

make_df <- function(res_obj) { as.data.frame(res_obj) %>% rownames_to_column("GeneID") }

write.csv(make_df(res_HET_vs_WT),  file.path(results_dir, "kallisto_gene_HET_vs_WT.csv"), row.names = FALSE)
write.csv(make_df(res_HOM_vs_WT),  file.path(results_dir, "kallisto_gene_HOM_vs_WT.csv"), row.names = FALSE)
write.csv(make_df(res_HOM_vs_HET), file.path(results_dir, "kallisto_gene_HOM_vs_HET.csv"), row.names = FALSE)

saveRDS(dds, file.path(results_dir, "dds_kallisto_gene.rds"))

