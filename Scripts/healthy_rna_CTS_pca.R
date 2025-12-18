## H2-0  Preparation ####

# Clean environment
rm(list = ls())

# Load required libraries
library(tidyverse)
library(sva)
library(edgeR)
library(ggplot2)
library(ggrepel)

## 1. Load your HEALTHY count data
cnt <- read.csv("Original_data/Combined_counts_Amsterdam_my_cohort_allhealthy.csv",
                check.names = FALSE)

# Remove Ensembl version numbers (e.g., ENSG00000123.5 → ENSG00000123)
ens_raw   <- as.character(cnt[[1]])
ens_clean <- sub("\\..*$", "", ens_raw)
cnt[[1]]  <- NULL
rownames(cnt) <- ens_clean

# Convert to numeric matrix
cnt_mat <- as.matrix(cnt)
storage.mode(cnt_mat) <- "numeric"

cat("Healthy cohort: Genes:", nrow(cnt_mat), " Samples:", ncol(cnt_mat), "\n")

## 2. Convert to logCPM
dge <- edgeR::DGEList(counts = cnt_mat)
expr_logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
cat("Healthy logCPM matrix dimensions:", nrow(expr_logcpm), "x", ncol(expr_logcpm), "\n")

## 3. Load Sepsis reference
load("Original_data/exp_core_g.rda")     # Expression of Sepsis reference samples
load("Original_data/core_samples.rda")   # Metadata of Sepsis reference

cat("Sepsis reference dimensions:", nrow(exp_core_g), "x", ncol(exp_core_g), "\n")

## 4. Create combined expression matrix (Sepsis ref + Healthy)
# Make sure genes overlap
common_genes <- intersect(rownames(exp_core_g), rownames(expr_logcpm))

expr_combined <- cbind(
  exp_core_g[common_genes, ],
  expr_logcpm[common_genes, ]
)

cat("Combined matrix dimensions (Sepsis + Healthy):",
    nrow(expr_combined), "x", ncol(expr_combined), "\n")

## 5. Create metadata automatically
meta <- data.frame(
  sample   = c(colnames(exp_core_g), colnames(expr_logcpm)),
  platform = c(rep("Sepsis_ref",  ncol(exp_core_g)),
               rep("Healthy_cohort", ncol(expr_logcpm))),
  cohort   = c(rep("Sepsis",        ncol(exp_core_g)),
               rep("Noninfectious", ncol(expr_logcpm)))
)

# Check metadata
head(meta)
table(meta$cohort)
table(meta$platform)

# Match sample order with expression columns
meta <- meta[match(colnames(expr_combined), meta$sample), ]
stopifnot(identical(colnames(expr_combined), meta$sample))

# Optional: save metadata
write.csv(meta, "Original_data/meta_platform_cohort_Healthy_vs_Sepsis.csv",
          row.names = FALSE)



## H2-1  Pre-ComBat visualization ####

# For density plot, randomly select up to 5000 genes
set.seed(1)
g_keep <- sample(rownames(expr_combined), min(5000, nrow(expr_combined)))
expr_sub <- expr_combined[g_keep, , drop = FALSE]

# --- Density plot (color = platform)
df_den <- as.data.frame(expr_sub) %>%
  dplyr::rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "logCPM") %>%
  dplyr::left_join(meta[, c("sample", "platform")], by = "sample")

p_den_pre <- ggplot(df_den, aes(x = logCPM, color = platform)) +
  geom_density() +
  theme_classic() +
  labs(
    title = "Pre-ComBat density by platform (Sepsis + Noninfectious)",
    x = "logCPM",
    y = "Density"
  )

# --- PCA (color = platform, shape = cohort)
pca_pre <- prcomp(t(expr_combined), center = TRUE, scale. = TRUE)

pca_df  <- data.frame(
  PC1 = pca_pre$x[,1],
  PC2 = pca_pre$x[,2],
  platform = meta$platform,
  cohort   = meta$cohort
)

p_pre <- ggplot(pca_df, aes(PC1, PC2, color = platform, shape = cohort)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_classic() +
  labs(
    title = "Pre-ComBat PCA (Sepsis + Noninfectious)",
    x = "PC1",
    y = "PC2"
  )

# Export
dir.create("Figure", showWarnings = FALSE)
ggsave("Figure/Pre_ComBat_PCA_noninfectious.svg", p_pre, width = 6, height = 5)
cat("✅ Saved: Figure/Pre_ComBat_PCA_noninfectious.svg\n")



## H2-2  Run ComBat batch correction ####

batch <- meta$platform
mod   <- model.matrix(~ cohort, data = meta)  # (not used below; kept for reference)

expr_combat <- ComBat(
  dat          = as.matrix(expr_combined),
  batch        = batch,
  mod          = NULL,   # consistent with your CAP script
  par.prior    = TRUE,
  prior.plots  = FALSE
)



## H2-3  Post-ComBat visualization ####

# --- Density plot after ComBat ---
expr_sub2 <- expr_combat[g_keep, , drop = FALSE]

df_den2 <- as.data.frame(expr_sub2) %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "logCPM") %>%
  left_join(meta[, c("sample","platform")], by = "sample")

p_den_post <- ggplot(df_den2, aes(x = logCPM, color = platform)) +
  geom_density() +
  theme_classic() +
  labs(
    title = "Post-ComBat density by platform (Sepsis + Noninfectious)",
    x = "logCPM",
    y = "Density"
  )

# --- PCA after ComBat ---
pca_post <- prcomp(t(expr_combat), center = TRUE, scale. = TRUE)

pca2_df  <- data.frame(
  PC1 = pca_post$x[,1],
  PC2 = pca_post$x[,2],
  platform = meta$platform,
  cohort   = meta$cohort
)

p_post <- ggplot(pca2_df, aes(PC1, PC2, color = platform, shape = cohort)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_classic() +
  labs(
    title = "Post-ComBat PCA (Sepsis + Noninfectious)",
    x = "PC1",
    y = "PC2"
  )

# Export
ggsave("Figure/Post_ComBat_PCA_noninfectious.svg",
       p_post, width = 6, height = 5)
cat("✅ Saved: Figure/Post_ComBat_PCA_noninfectious.svg\n")
