rm(list = ls())

# ---- Packages ----
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
cran_pkgs <- c("dplyr","readr","tidyr","stringr","ggplot2","RColorBrewer","tableone")
for(p in cran_pkgs) if(!requireNamespace(p, quietly = TRUE)) install.packages(p)
bio_pkgs <- c("biomaRt")
for(p in bio_pkgs) if(!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE)
library(dplyr); library(readr); library(tidyr); library(stringr)
library(ggplot2); library(RColorBrewer); library(tableone); library(biomaRt)
library(readr)
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
               tableone, kmed, survival, survminer, skimr, lattice, reshape2, grDevices, ggrepel)

# ready data
# use the "Mixture_for_CIBERSORTx_RNAseq_TPM.csv" made by code MS1
# then go to the website CIBERSORTx
# setting of the website of my data: Job type: Impute Cell Fractions; Signature matrix file: K999 signature (filename above)
   #Batch correction: enabled, S-mode; S-mode reference: ref_matrix_equal.tsv; Quantile normalization: Disable quantile normalization (disabling is recommended for RNA-Seq data)
   #Run mode: absolute; Permutations: 100

# make the percentage for IL1R2
# ===============================
# Calculate IL1R2 signature (absolute and percentage)
# ===============================

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)

# --- Step 1. Read full CIBERSORTx output ---
cib_path <- "Original_data/CIBERSORTx_Job3_Adjusted_IL1R2.txt"
x <- read.delim(cib_path, check.names = FALSE)
colnames(x)[colnames(x) == "IL1R2+_immature_neutrophils"] <- "IL1R2"

# --- Step 2. Define metadata and cell-state columns ---
meta <- c("Mixture", "P-value", "Correlation", "RMSE", "Absolute score (sig.score)")
state_cols <- setdiff(names(x), meta)

# ===============================
# A) IL1R2 absolute abundance
# ===============================
# absolute abundance = fraction * Absolute score
il1r2_abs <- x %>%
  transmute(
    Mixture,
    IL1R2_absolute = IL1R2 * `Absolute score (sig.score)`
  )

out_abs <- "Original_data/CIBERSORTx_IL1R2_absolute.csv"
write_csv(il1r2_abs, out_abs)
message("Saved absolute IL1R2 abundance: ", out_abs)

# ===============================
# B) IL1R2 percentage (within-sample composition)
# ===============================
# Use sum across all cell-state columns as the total, then compute IL1R2%
row_tot <- rowSums(x[, state_cols, drop = FALSE], na.rm = TRUE)

il1r2_pct <- x %>%
  transmute(
    Mixture,
    IL1R2_percent = IL1R2 / row_tot * 100
  )

out_pct <- "Original_data/CIBERSORTx_IL1R2_percent.csv"
write_csv(il1r2_pct, out_pct)
message("Saved IL1R2 percentage: ", out_pct)

# --- QC check: each row should sum ≈ 100 ---
pct_all <- sweep(x[, state_cols, drop = FALSE], 1, row_tot, "/") * 100
cat("QC (row-wise % sum - 100) summary:\n")
print(summary(rowSums(pct_all) - 100))

# ===============================
# C) Merge both for inspection
# ===============================
out_both <- il1r2_abs %>%
  left_join(il1r2_pct, by = "Mixture") %>%
  mutate(Absolute_score = x$`Absolute score (sig.score)`)

out_both_path <- "Original_data/CIBERSORTx_IL1R2_absolute_and_percent.csv"
write_csv(out_both, out_both_path)
message("Saved: ", out_both_path)

# ===============================
# D) Visualization
# ===============================

# Histogram: absolute abundance
p_abs <- ggplot(il1r2_abs, aes(x = IL1R2_absolute)) +
  geom_histogram(bins = 30, fill = "#E69F00", color = "white", alpha = 0.8) +
  geom_vline(aes(xintercept = median(IL1R2_absolute, na.rm = TRUE)),
             color = "red", linetype = "dashed", size = 1) +
  labs(title = "IL1R2 absolute abundance (Absolute mode, unnormalized)",
       x = "IL1R2 absolute score", y = "Count") +
  theme_classic(base_size = 14)
print(p_abs)

# Histogram: percentage
p_pct <- ggplot(il1r2_pct, aes(x = IL1R2_percent)) +
  geom_histogram(bins = 30, fill = "#56B4E9", color = "white", alpha = 0.8) +
  geom_vline(aes(xintercept = median(IL1R2_percent, na.rm = TRUE)),
             color = "red", linetype = "dashed", size = 1) +
  labs(title = "IL1R2 percentage (within-sample composition)",
       x = "IL1R2 proportion (%)", y = "Count") +
  theme_classic(base_size = 14)
print(p_pct)
