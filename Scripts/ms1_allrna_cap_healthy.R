## ============================================================
## 0) CLEAN
## ============================================================
rm(list = ls())

## ============================================================
## 1) Load packages
## ============================================================
if(!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
if(!requireNamespace("biomaRt", quietly=TRUE)) BiocManager::install("biomaRt")

library(dplyr)
library(readr)
library(biomaRt)

## ============================================================
## 2) Read two raw count matrices (Ensembl IDs in first column)
## ============================================================
cnt_healthy <- read.csv("Original_data/Combined_counts_Amsterdam_my_cohort_allhealthy.csv",
                        check.names = FALSE)

cnt_cap     <- read.csv("Original_data/Combined_counts_Amsterdam_my_cohort_allpatients.csv",
                        check.names = FALSE)

## Check dimensions
cat("Healthy:",  dim(cnt_healthy), "\n")
cat("CAP:",      dim(cnt_cap), "\n")

## ============================================================
## 3) Combine TWO datasets into one big matrix (safe version)
## ============================================================

# put Ensembl IDs as rownames
ens_healthy <- sub("\\..*$", "", cnt_healthy[[1]])
cnt_healthy[[1]] <- NULL
rownames(cnt_healthy) <- ens_healthy

ens_cap <- sub("\\..*$", "", cnt_cap[[1]])
cnt_cap[[1]] <- NULL
rownames(cnt_cap) <- ens_cap

# make sure all columns are numeric
cnt_healthy <- as.data.frame(lapply(cnt_healthy, as.numeric))
cnt_cap     <- as.data.frame(lapply(cnt_cap, as.numeric))
rownames(cnt_healthy) <- ens_healthy
rownames(cnt_cap)     <- ens_cap

# get all unique genes
all_genes <- union(rownames(cnt_healthy), rownames(cnt_cap))

# full matrix initialized with zeros
full_cnt <- matrix(0,
                   nrow = length(all_genes),
                   ncol = ncol(cnt_healthy) + ncol(cnt_cap))
rownames(full_cnt) <- all_genes
colnames(full_cnt) <- c(colnames(cnt_healthy), colnames(cnt_cap))

# fill in
full_cnt[rownames(cnt_healthy), colnames(cnt_healthy)] <- as.matrix(cnt_healthy)
full_cnt[rownames(cnt_cap),     colnames(cnt_cap)]     <- as.matrix(cnt_cap)

cat("Full combined matrix:", nrow(full_cnt), "genes ×", ncol(full_cnt), "samples\n")

# convert to data.frame for downstream TPM calculation
full_cnt <- as.data.frame(full_cnt)

## ============================================================
## 4) Fetch gene lengths for Ensembl IDs
## ============================================================
mart <- useEnsembl(biomart="genes", dataset="hsapiens_gene_ensembl")

annot <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "start_position", "end_position"),
  filters    = "ensembl_gene_id",
  values     = rownames(full_cnt),
  mart       = mart
)

annot <- annot %>%
  mutate(gene_length = end_position - start_position + 1) %>%
  filter(!is.na(gene_length) & gene_length > 0) %>%
  group_by(ensembl_gene_id) %>%
  slice(1) %>% ungroup()

## align
common <- intersect(rownames(full_cnt), annot$ensembl_gene_id)
cnt_use <- full_cnt[common, ]
ann_use <- annot[match(common, annot$ensembl_gene_id), ]

## ============================================================
## 5) Compute TPM
## ============================================================
cat("Computing TPM...\n")

length_kb <- ann_use$gene_length / 1000
rate      <- sweep(cnt_use, 1, length_kb, "/")   # counts per kb

# Collapse to HGNC symbol
## ============================================================
## Collapse Ensembl → HGNC symbol (sum rates)
## ============================================================

rate_df <- cbind(symbol = ann_use$hgnc_symbol, as.data.frame(rate))

# remove NA/empty symbol
rate_df <- rate_df[!is.na(rate_df$symbol) & rate_df$symbol != "", ]

# split by symbol
spl <- split(rate_df[,-1], rate_df$symbol)

# sum over Ensembl rows for each symbol
rate2_mat <- sapply(spl, function(x) colSums(x, na.rm = TRUE))

# transpose to symbol × sample
rate2_mat <- t(rate2_mat)

# final symbol list
symbols <- rownames(rate2_mat)

cat("Collapsed to", length(symbols), "unique HGNC symbols\n")

rate_mat <- rate2_mat  # rename to keep consistent with your code

# scale to TPM
scaler  <- colSums(rate_mat) / 1e6
tpm_mat <- sweep(rate_mat, 2, scaler, "/")

tpm_df <- cbind(symbol = symbols, as.data.frame(tpm_mat))
tpm_df <- tpm_df[rowSums(tpm_df[,-1]) > 0, ]

cat("TPM ready. Genes:", nrow(tpm_df), "\n")

## ============================================================
## 6) Export TWO sets:
##    (1) TPMS for full cohort (healthy + CAP)
##    (2) TPM only for healthy (non-infection)
## ============================================================

# get sample names
all_samples <- colnames(tpm_df)[-1]
healthy_samples <- colnames(cnt_healthy)
cap_samples     <- colnames(cnt_cap)

# (1) TPM for FULL cohort
full_tpm <- tpm_df[, c("symbol", all_samples)]
write_csv(full_tpm, "Original_data/TPM_full_noninf_CAP.csv")
readr::write_tsv(full_tpm, "Original_data/TPM_full_noninf_CAP.txt", na = "")

# (2) TPM only for healthy
healthy_intersect <- intersect(all_samples, healthy_samples)
healthy_tpm <- tpm_df[, c("symbol", healthy_intersect)]
write_csv(healthy_tpm, "Original_data/TPM_only_noninfection.csv")
readr::write_tsv(healthy_tpm, "Original_data/TPM_only_noninfection.txt", na = "")




#after CIBERSORTx
## ============================================================
## Load packages
## ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

## ============================================================
## 1. Read CIBERSORTx results
## ============================================================

file_non     <- "Original_data/CIBERSORTx_non_Results.txt"
file_non_cap <- "Original_data/CIBERSORTx_non_cap_Results.txt"

cib_non     <- read.delim(file_non,     check.names = FALSE)
cib_non_cap <- read.delim(file_non_cap, check.names = FALSE)

## ============================================================
## 2. Extract correct MS1 values (NO multiplying with sig.score)
## ============================================================

ms1_non <- cib_non %>%
  mutate(MS1_non = MS1) %>%
  select(Mixture, MS1_non)

ms1_noncap <- cib_non_cap %>%
  mutate(MS1_noncap = MS1) %>%
  select(Mixture, MS1_noncap)

## ============================================================
## 3. Merge two MS1 sets
## ============================================================

ms1_merged <- ms1_non %>%
  inner_join(ms1_noncap, by = "Mixture")

## ============================================================
## 4. Scatter plot
## ============================================================

cor_value <- round(cor(ms1_merged$MS1_non,
                       ms1_merged$MS1_noncap,
                       use = "complete.obs"), 3)

p_scatter <- ggplot(ms1_merged,
                    aes(x = MS1_non, y = MS1_noncap)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 1) +
  geom_smooth(method = "lm", se = FALSE, color = "blue", linewidth = 0.8) +
  annotate("text",
           x = min(ms1_merged$MS1_non, na.rm = TRUE),
           y = max(ms1_merged$MS1_noncap, na.rm = TRUE),
           hjust = 0,
           label = paste("Pearson r =", cor_value),
           size = 5) +
  labs(
    title = "MS1 for non-infection samples under different model inputs",
    x = "MS1 (NON-only input)",
    y = "MS1 (NON subset from NON+CAP input)"
  ) +
  theme_classic(base_size = 16)

print(p_scatter)

## ============================================================
## 5. Density plot
## ============================================================

abs_long <- ms1_merged %>%
  pivot_longer(cols = c(MS1_non, MS1_noncap),
               names_to = "InputType",
               values_to = "MS1") %>%
  mutate(InputType = factor(InputType,
                            levels = c("MS1_non", "MS1_noncap"),
                            labels = c("NON-only input",
                                       "NON subset from NON+CAP input")))

p_density_line <- ggplot(abs_long,
                         aes(x = MS1,
                             color = InputType,
                             linetype = InputType)) +
  geom_density(linewidth = 1.2) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  scale_color_manual(values = c("#1b9e77", "#d95f02")) +
  labs(
    title = "MS1 distribution of non-infection samples",
    x = "MS1",
    y = "Density"
  ) +
  theme_classic(base_size = 16) +
  theme(
    legend.title = element_blank(),
    legend.position = "top"
  )

print(p_density_line)

## ============================================================
## 6. Boxplot
## ============================================================

p_box <- ggplot(abs_long,
                aes(x = InputType, y = MS1, fill = InputType)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = c("#1b9e77", "#d95f02")) +
  labs(
    title = "MS1 values for non-infection samples",
    x = "",
    y = "MS1"
  ) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12)
  )

print(p_box)


## ================================
## Load packages
## ================================
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

## ================================
## 1. Read CIBERSORTx results
## ================================
file_non     <- "Original_data/CIBERSORTx_non_Results.txt"      # only NON
file_non_cap <- "Original_data/CIBERSORTx_non_cap_Results.txt"  # NON + CAP

cib_non     <- read.delim(file_non,     check.names = FALSE)
cib_non_cap <- read.delim(file_non_cap, check.names = FALSE)

## ================================
## 2. 定义 NON 样本的 ID（来自 non-only 那个文件）
## ================================
non_ids <- unique(cib_non$Mixture)

## ================================
## 3. 在 non+CAP 这次 run 里面计算 MS1 percentage
##    MS1% = MS1 / sum(16 signatures) * 100
## ================================
meta_cols_cap <- c("Mixture",
                   "P-value", "P.value",          # 可能的名字二选一
                   "Correlation",
                   "RMSE",
                   "Absolute score (sig.score)",
                   "Absolute.score..sig.score.")

state_cols_cap <- setdiff(names(cib_non_cap), meta_cols_cap)

cib_pct <- cib_non_cap %>%
  mutate(
    row_total    = rowSums(across(all_of(state_cols_cap)), na.rm = TRUE),
    MS1_percent  = MS1 / row_total * 100,
    Group = ifelse(Mixture %in% non_ids, "NON-infection", "CAP")
  ) %>%
  select(Mixture, Group, MS1_percent)

table(cib_pct$Group)   # 检查一下 CAP 和 NON-infection 是否都出现

## ================================
## 4. Density plot: CAP vs NON-infection
## ================================
p_density <- ggplot(cib_pct,
                    aes(x = MS1_percent, color = Group)) +
  geom_density(size = 1.2) +
  scale_color_manual(values = c("NON-infection" = "#1b9e77",
                                "CAP"           = "#d95f02")) +
  labs(title = "MS1 percentage distribution: CAP vs non-infection",
       x = "MS1 percentage (%)",
       y = "Density") +
  theme_classic(base_size = 16) +
  theme(legend.title = element_blank(),
        legend.position = "top")

print(p_density)
# ggsave("MS1_pct_density_CAP_vs_NON.pdf", p_density, width = 7, height = 5)

## ================================
## 5. Boxplot: CAP vs NON-infection
## ================================
p_box <- ggplot(cib_pct,
                aes(x = Group, y = MS1_percent, fill = Group)) +
  geom_boxplot(alpha = 0.8, width = 0.6) +
  scale_fill_manual(values = c("NON-infection" = "#1b9e77",
                               "CAP"           = "#d95f02")) +
  labs(title = "MS1 percentage: CAP vs non-infection",
       x = "",
       y = "MS1 percentage (%)") +
  theme_classic(base_size = 16) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12))

print(p_box)
# ggsave("MS1_pct_boxplot_CAP_vs_NON.pdf", p_box, width = 6, height = 5)




#percentage
## ================================
## Load packages
## ================================
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

## ================================
## 1. Read CIBERSORTx results
## ================================
file_non     <- "Original_data/CIBERSORTx_non_Results.txt"      # only NON
file_non_cap <- "Original_data/CIBERSORTx_non_cap_Results.txt"  # NON + CAP

cib_non     <- read.delim(file_non,     check.names = FALSE)
cib_non_cap <- read.delim(file_non_cap, check.names = FALSE)

## ================================
## 2. 定义 NON 样本的 ID（来自 non-only 那个文件）
## ================================
non_ids <- unique(cib_non$Mixture)

## ================================
## 3. 在 non+CAP 这次 run 里面计算 MS1 percentage
##    MS1% = MS1 / sum(16 signatures) * 100
## ================================
meta_cols_cap <- c("Mixture",
                   "P-value", "P.value",          # 可能的名字二选一
                   "Correlation",
                   "RMSE",
                   "Absolute score (sig.score)",
                   "Absolute.score..sig.score.")

state_cols_cap <- setdiff(names(cib_non_cap), meta_cols_cap)

cib_pct <- cib_non_cap %>%
  mutate(
    row_total    = rowSums(across(all_of(state_cols_cap)), na.rm = TRUE),
    MS1_percent  = MS1 / row_total * 100,
    Group = ifelse(Mixture %in% non_ids, "NON-infection", "CAP")
  ) %>%
  select(Mixture, Group, MS1_percent)

table(cib_pct$Group)   # 检查一下 CAP 和 NON-infection 是否都出现

## ================================
## 4. Density plot: CAP vs NON-infection
## ================================
p_density <- ggplot(cib_pct,
                    aes(x = MS1_percent, color = Group)) +
  geom_density(size = 1.2) +
  scale_color_manual(values = c("NON-infection" = "#1b9e77",
                                "CAP"           = "#d95f02")) +
  labs(title = "MS1 signature proportion: CAP vs non-infection",
       x = "MS1 signature proportion (%)",
       y = "Density") +
  theme_classic(base_size = 16) +
  theme(legend.title = element_blank(),
        legend.position = "top")

print(p_density)
# ggsave("MS1_pct_density_CAP_vs_NON.pdf", p_density, width = 7, height = 5)

## ================================
## 5. Boxplot: CAP vs NON-infection
## ================================
p_box <- ggplot(cib_pct,
                aes(x = Group, y = MS1_percent, fill = Group)) +
  geom_boxplot(alpha = 0.8, width = 0.6) +
  scale_fill_manual(values = c("NON-infection" = "#1b9e77",
                               "CAP"           = "#d95f02")) +
  labs(title = "MS1 signature proportion: CAP vs non-infection",
       x = "",
       y = "MS1 percentage (%)") +
  theme_classic(base_size = 16) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12))

print(p_box)
# ggsave("MS1_pct_boxplot_CAP_vs_NON.pdf", p_box, width = 6, height = 5)
