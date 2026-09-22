# High-resolution transcript-level differential expression analysis using Sleuth

# Script 07: Uncollapsed Transcript-Level Statistics via Sleuth
# Statistical Framework: Response Models Accounting for Inferential Variance

library(devtools, lib.loc="~/Rlib")
library(sleuth, lib.loc="~/Rlib")
library(dplyr)

base_dir <- "/users/y3/ag2622/kallisto_tximport"
out_dir  <- "/users/y3/ag2622/new/kallisto_sleuth"

sample_id <- c("SRR8985047", "SRR8985048", "SRR8985049", "SRR8985050", "SRR8985051", "SRR8985052")
condition <- c("WT", "WT", "HET", "HET", "HOM", "HOM")
kal_dirs  <- file.path(base_dir, "quant", sample_id)

s2c <- data.frame(sample = sample_id, condition = condition, path = kal_dirs, stringsAsFactors = FALSE)
s2c$condition <- factor(s2c$condition, levels = c("WT", "HET", "HOM"))

# Prep loads the 100 internal bootstrap variations
so <- sleuth_prep(s2c, ~ condition)
so <- sleuth_fit(so, ~ condition, 'full')
so <- sleuth_fit(so, ~ 1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')

# Wald tests to bypass label compilation errors
so <- sleuth_wt(so, 'conditionHET')
so <- sleuth_wt(so, 'conditionHOM')

res_HET_vs_WT <- sleuth_results(so, 'conditionHET')
res_HOM_vs_WT <- sleuth_results(so, 'conditionHOM')

# Extract HET reference contrast shift matrix
s2c_relevel <- s2c
s2c_relevel$condition <- relevel(s2c_relevel$condition, "HET")
so2            <- sleuth_prep(s2c_relevel, ~ condition)
so2            <- sleuth_fit(so2, ~ condition, 'full')
so2            <- sleuth_wt(so2, 'conditionHOM')
res_HOM_vs_HET <- sleuth_results(so2, 'conditionHOM')

filter_sleuth <- function(df) { df %>% filter(!is.na(qval) & qval < 0.05 & abs(b) > 1) }

write.csv(as.data.frame(filter_sleuth(res_HET_vs_WT)),  file.path(out_dir, "sleuth_HET_vs_WT_sig.csv"), row.names = FALSE)
write.csv(as.data.frame(filter_sleuth(res_HOM_vs_WT)),  file.path(out_dir, "sleuth_HOM_vs_WT_sig.csv"), row.names = FALSE)
write.csv(as.data.frame(filter_sleuth(res_HOM_vs_HET)), file.path(out_dir, "sleuth_HOM_vs_HET_sig.csv"), row.names = FALSE)

