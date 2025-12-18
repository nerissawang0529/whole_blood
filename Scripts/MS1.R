# ===============================
# MS1 with CIBERSORTx: full recipe
# From raw counts (Ensembl) -> TPM (HGNC) -> CIBERSORTx -> MS1 -> Table 1
# ===============================

rm(list = ls())

# ---- Packages ----
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
cran_pkgs <- c("dplyr","readr","tidyr","stringr","ggplot2","RColorBrewer","tableone")
for(p in cran_pkgs) if(!requireNamespace(p, quietly = TRUE)) install.packages(p)
bio_pkgs <- c("biomaRt")
for(p in bio_pkgs) if(!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE)

library(dplyr); library(readr); library(tidyr); library(stringr)
library(ggplot2); library(RColorBrewer); library(tableone); library(biomaRt)

# -------------------------------------------------------------------
# 1) READ COUNTS (Ensembl IDs in column 1, samples in columns 2..N)
# -------------------------------------------------------------------
counts_file <- "Original_data/Combined_counts_Amsterdam_my_cohort.csv"
cnt <- read.csv(counts_file, check.names = FALSE)
stopifnot(ncol(cnt) > 2)

# strip Ensembl version (e.g., ENSG... .12 -> ENSG...)
ens_raw   <- as.character(cnt[[1]])
ens_novers <- sub("\\..*$", "", ens_raw)
cnt[[1]]  <- NULL
rownames(cnt) <- ens_novers

# numeric matrix
cnt_mat <- as.matrix(cnt)
storage.mode(cnt_mat) <- "numeric"

# collapse duplicate Ensembl IDs (after version stripping) by summing counts
if (any(duplicated(rownames(cnt_mat)))) {
  cnt_mat <- rowsum(cnt_mat, group = rownames(cnt_mat))  # sum rows with same name
}

cat("Counts matrix -> Genes:", nrow(cnt_mat), " Samples:", ncol(cnt_mat), "\n")

# -------------------------------------------------------------------
# 2) GET GENE LENGTHS (bp) + MAP ENSEMBL -> HGNC SYMBOLS
# -------------------------------------------------------------------
message("Fetching gene coords & HGNC symbols from Ensembl (biomaRt)...")

# You can choose a mirror if needed: mirror = "www" | "uswest" | "useast" | "asia"
mart <- biomaRt::useEnsembl(biomart = "genes",
                            dataset = "hsapiens_gene_ensembl")

attrs <- c("ensembl_gene_id", "hgnc_symbol", "start_position", "end_position")
annot_raw <- biomaRt::getBM(
  attributes = attrs,
  filters    = "ensembl_gene_id",
  values     = rownames(cnt_mat),
  mart       = mart
)

# Compute gene length (bp) from genomic coords
annot <- annot_raw |>
  dplyr::mutate(gene_length = as.numeric(end_position) - as.numeric(start_position) + 1L) |>
  dplyr::select(ensembl_gene_id, hgnc_symbol, gene_length) |>
  dplyr::filter(!is.na(gene_length), gene_length > 0)

# Keep one row per Ensembl (prefer non-empty HGNC)
annot <- annot |>
  dplyr::group_by(ensembl_gene_id) |>
  dplyr::arrange(dplyr::desc(nchar(hgnc_symbol))) |>
  dplyr::slice(1) |>
  dplyr::ungroup()


# align & merge
common_ens <- intersect(rownames(cnt_mat), annot$ensembl_gene_id)
stopifnot(length(common_ens) > 500)  # sanity check

cnt_use <- cnt_mat[common_ens, , drop = FALSE]
ann_use <- annot[match(common_ens, annot$ensembl_gene_id), ]

# collapse to unique HGNC symbols: sum counts, median length
dat <- cbind(hgnc_symbol = ann_use$hgnc_symbol,
             gene_length = ann_use$gene_length,
             as.data.frame(cnt_use, check.names = FALSE))

# ===== Correct collapse -> TPM =====
# cnt_use: rows = Ensembl IDs, cols = samples
# ann_use: columns: ensembl_gene_id, hgnc_symbol, gene_length (bp)

# 1) counts per kilobase for each Ensembl row
length_kb <- ann_use$gene_length / 1000
rate <- sweep(cnt_use, 1, length_kb, "/")   # rate = counts / kb

# 2) attach SYMBOL and collapse by SYMBOL by summing rates
if ("package:plyr" %in% search()) detach("package:plyr", unload = TRUE)
library(dplyr)

rate_dat <- cbind(hgnc_symbol = ann_use$hgnc_symbol,
                  as.data.frame(rate, check.names = FALSE)) |>
  tibble::as_tibble() |>
  filter(!is.na(hgnc_symbol) & hgnc_symbol != "") |>
  group_by(hgnc_symbol) |>
  summarise(
    across(.cols = where(is.numeric), .fns = ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  )


# 3) scale to TPM per sample
rate_mat <- as.matrix(rate_dat[,-1, drop = FALSE])         # rows = SYMBOL, cols = samples
scaler   <- colSums(rate_mat, na.rm = TRUE) / 1e6          # sum(rate) / 1e6
tpm_sym  <- sweep(rate_mat, 2, scaler, "/")                # TPM

# QC: TPM columns should sum to ~1e6
tpmsums <- colSums(tpm_sym)
cat("Median TPM column sum:", median(tpmsums), "\n")

# 4) final Mixture for CIBERSORTx (RNA-seq TPM, non-log)
tpm_df <- cbind(
  data.frame(gennames = rate_dat$hgnc_symbol, check.names = FALSE),
  as.data.frame(tpm_sym, check.names = FALSE)
)

# (optional) drop all-zero rows
tpm_df <- tpm_df[rowSums(tpm_df[,-1, drop = FALSE]) > 0, ]

# write file
out_csv <- file.path("Original_data", "Mixture_for_CIBERSORTx_RNAseq_TPM.csv")
readr::write_csv(tpm_df, out_csv)
message("Saved: ", out_csv)


out_tsv <- file.path("Original_data","Mixture_for_CIBERSORTx_RNAseq_TPM.txt")
readr::write_tsv(tpm_df, out_tsv, na = "")   # readr never adds quotes; NA written as empty (but we asserted none)
message("Saved: ", out_tsv)

##then go to the website CIBERSORTx####

##the get the MS1 score

## Calculate percentage of each Cell_stage for each patient and Stacked barplot
# Giuseppe Leite, 
# giuseppe.gianini@unifesp.br
# Amsterdam UMC, locatie AMC, April 2023
# Center of Experimental & Molecular Medicine (CEMM)
# setting of the website of my data: 2.impute cell fractions, Custom, Enable batch correction (B-mode),Disable quantile normalization, Run in absolute mode, Permutations for significance analysis 100 


# ===============================
# Calculate absolute abundance for each cell state (no within-sample normalization)
# ===============================

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)

# Read the full CIBERSORTx output (do NOT subset columns like select(1:17))
cib_path <- "Original_data/CIBERSORTx_Job11_Adjusted.txt"
CIBERSORTx_results <- read.delim(cib_path, check.names = FALSE)

# Define metadata columns (adapted to your current file header)
meta <- c("Mixture","P-value","Correlation","RMSE","Absolute score (sig.score)")
state_cols <- setdiff(names(CIBERSORTx_results), meta)  # all deconvolved cell-state columns

# Compute per–cell-state absolute abundance = fraction * Absolute score
CIBERSORTx_absolute <- CIBERSORTx_results %>%
  dplyr::mutate(across(all_of(state_cols), ~ .x * `Absolute score (sig.score)`)) %>%
  dplyr::select(Mixture, all_of(state_cols))

# QC: row sums of absolute matrix should ≈ Absolute score
check_sum <- rowSums(CIBERSORTx_absolute[,-1, drop = FALSE])
qc_diff   <- check_sum - CIBERSORTx_results$`Absolute score (sig.score)`
cat("QC (absolute row sum - Absolute score) summary:\n")
print(summary(qc_diff))

# Extract MS1 absolute abundance and save
CIBERSORTx_absolute_MS1 <- CIBERSORTx_absolute %>%
  dplyr::select(Mixture, MS1)

out_abs_ms1 <- file.path("Original_data", "CIBERSORTx_absolute_MS1.csv")
write_csv(CIBERSORTx_absolute_MS1, out_abs_ms1)
message("Saved absolute MS1 abundance: ", out_abs_ms1)

# Plot histogram of MS1 absolute abundance
p_abs <- ggplot(CIBERSORTx_absolute_MS1, aes(x = MS1)) +
  geom_histogram(bins = 30, fill = "#E69F00", color = "white", alpha = 0.8) +
  geom_vline(aes(xintercept = median(MS1, na.rm = TRUE)),
             color = "red", linetype = "dashed", size = 1) +
  labs(title = "MS1 absolute abundance (Absolute mode, unnormalized)",
       x = "MS1 absolute score",
       y = "Count") +
  theme_classic(base_size = 14)
print(p_abs)
# ggsave("Original_data/MS1_absolute_histogram.pdf", p_abs, width = 6, height = 4)


# ===============================
# Derive both absolute and percentage metrics from the same file
# ===============================

# Re-read (or reuse) the full output to avoid accidental column truncation
x <- read.delim(cib_path, check.names = FALSE)

# Metadata and state columns (same as above)
meta <- c("Mixture","P-value","Correlation","RMSE","Absolute score (sig.score)")
state_cols <- setdiff(names(x), meta)

## ---------- A) MS1 absolute abundance ----------
ms1_abs <- x %>%
  transmute(
    Mixture,
    MS1_absolute = MS1 * `Absolute score (sig.score)`
  )
out_abs <- "Original_data/CIBERSORTx_MS1_absolute.csv"
write_csv(ms1_abs, out_abs)
message("Saved: ", out_abs)

## ---------- B) MS1 percentage (within-sample composition) ----------
# Use the sum across all state columns as the row total, then compute MS1%
row_tot <- rowSums(x[, state_cols, drop = FALSE], na.rm = TRUE)
ms1_pct <- x %>%
  transmute(
    Mixture,
    MS1_percent = MS1 / row_tot * 100
  )
out_pct <- "Original_data/CIBERSORTx_MS1_percent.csv"
write_csv(ms1_pct, out_pct)
message("Saved: ", out_pct)

# Optional: quick QC that percentages across states sum to ~100 for each sample
pct_all <- sweep(x[, state_cols, drop = FALSE], 1, row_tot, "/") * 100
cat("QC (row-wise % sum - 100) summary:\n")
print(summary(rowSums(pct_all) - 100))

# (Optional) combine both metrics for side-by-side inspection
out_both <- ms1_abs %>%
  left_join(ms1_pct, by = "Mixture") %>%
  mutate(Absolute_score = x$`Absolute score (sig.score)`)
out_both_path <- "Original_data/CIBERSORTx_MS1_absolute_and_percent.csv"
write_csv(out_both, out_both_path)
message("Saved: ", out_both_path)

# Optional: plot MS1% distribution
p_pct <- ggplot(ms1_pct, aes(x = MS1_percent)) +
  geom_histogram(bins = 30, fill = "#56B4E9", color = "white", alpha = 0.8) +
  geom_vline(aes(xintercept = median(MS1_percent, na.rm = TRUE)),
             color = "red", linetype = "dashed", size = 1) +
  labs(title = "MS1 percentage (within-sample composition)",
       x = "MS1 proportion (%)",
       y = "Count") +
  theme_classic(base_size = 14)
print(p_pct)
# ggsave("Original_data/MS1_percent_histogram.pdf", p_pct, width = 6, height = 4)
