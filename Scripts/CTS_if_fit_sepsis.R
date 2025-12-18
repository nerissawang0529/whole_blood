## H1-0  Preparation ####

# Clean environment
rm(list = ls())

# Load required libraries
library(tidyverse)
library(sva)
library(edgeR)
library(ggplot2)
library(ggrepel)

#1. Load your CAP count data
cnt <- read.csv("Original_data/Combined_counts_Amsterdam_my_cohort_allpatients.csv",
                check.names = FALSE)

# Remove Ensembl version numbers (e.g., ENSG00000123.5 → ENSG00000123)
ens_raw   <- as.character(cnt[[1]])
ens_clean <- sub("\\..*$", "", ens_raw)
cnt[[1]]  <- NULL
rownames(cnt) <- ens_clean

# Convert to numeric matrix
cnt_mat <- as.matrix(cnt)
storage.mode(cnt_mat) <- "numeric"

cat("Genes:", nrow(cnt_mat), " Samples:", ncol(cnt_mat), "\n")

## 2. Convert to logCPM
dge <- edgeR::DGEList(counts = cnt_mat)
expr_logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
cat("logCPM matrix dimensions:", nrow(expr_logcpm), "x", ncol(expr_logcpm), "\n")

## 3. Load Sepsis reference
load("Original_data/exp_core_g.rda")     # Expression of Sepsis reference samples
load("Original_data/core_samples.rda")   # Metadata of Sepsis reference

cat("Sepsis reference dimensions:", nrow(exp_core_g), "x", ncol(exp_core_g), "\n")

## 4. Create combined expression matrix
# Make sure genes overlap
common_genes <- intersect(rownames(exp_core_g), rownames(expr_logcpm))
expr_combined <- cbind(exp_core_g[common_genes, ],
                       expr_logcpm[common_genes, ])

cat("Combined matrix dimensions:", nrow(expr_combined), "x", ncol(expr_combined), "\n")

##5. Create metadata automatically
meta <- data.frame(
  sample   = c(colnames(exp_core_g), colnames(expr_logcpm)),
  platform = c(rep("Sepsis_ref", ncol(exp_core_g)),
               rep("CAP_cohort", ncol(expr_logcpm))),
  cohort   = c(rep("Sepsis", ncol(exp_core_g)),
               rep("CAP", ncol(expr_logcpm)))
)

# Check metadata
head(meta)
table(meta$cohort)
table(meta$platform)

# Match sample order with expression columns
meta <- meta[match(colnames(expr_combined), meta$sample), ]
stopifnot(identical(colnames(expr_combined), meta$sample))

# Optional: save metadata
write.csv(meta, "Original_data/meta_platform_cohort.csv", row.names = FALSE)



## H1-1  Pre-ComBat visualization


# For density plot, randomly select up to 5000 genes
set.seed(1)
g_keep <- sample(rownames(expr_combined), min(5000, nrow(expr_combined)))
expr_sub <- expr_combined[g_keep, , drop = FALSE]

# --- Density plot (color = platform)
df_den <- as.data.frame(expr_sub) %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "logCPM") %>%
  left_join(meta[, c("sample", "platform")], by = "sample")

ggplot(df_den, aes(x = logCPM, color = platform)) +
  geom_density() +
  theme_classic() +
  labs(title = "Pre-ComBat density by platform", x = "logCPM", y = "Density")

# --- PCA (color = platform, shape = cohort)
pca_pre <- prcomp(t(expr_combined), center = TRUE, scale. = TRUE)
pca_df  <- data.frame(PC1 = pca_pre$x[,1],
                      PC2 = pca_pre$x[,2],
                      platform = meta$platform,
                      cohort = meta$cohort)

p_pre <- ggplot(pca_df, aes(PC1, PC2, color = platform, shape = cohort)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_classic() +
  labs(title = "Pre-ComBat PCA", x = "PC1", y = "PC2")

# Export
dir.create("Figure", showWarnings = FALSE)
ggsave("Figure/Pre_ComBat_PCA.svg", p_pre, width = 6, height = 5)
cat("✅ Saved: Figure/Pre_ComBat_PCA.svg\n")



## H1-2  Run ComBat batch correction


batch <- meta$platform
mod   <- model.matrix(~ cohort, data = meta)  # Keep CAP vs Sepsis difference

expr_combat <- ComBat(
  dat          = as.matrix(expr_combined),
  batch        = batch,
  mod          = NULL,   # <-- remove this line
  par.prior    = TRUE,
  prior.plots  = FALSE
)




## H1-3  Post-ComBat visualization


# --- Density plot after ComBat ---
expr_sub2 <- expr_combat[g_keep, , drop = FALSE]

df_den2 <- as.data.frame(expr_sub2) %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "logCPM") %>%
  left_join(meta[, c("sample","platform")], by = "sample")

ggplot(df_den2, aes(x = logCPM, color = platform)) +
  geom_density() +
  theme_classic() +
  labs(title = "Post-ComBat density by platform", x = "logCPM", y = "Density")

# --- PCA after ComBat ---
pca_post <- prcomp(t(expr_combat), center = TRUE, scale. = TRUE)
pca2_df  <- data.frame(PC1 = pca_post$x[,1],
                       PC2 = pca_post$x[,2],
                       platform = meta$platform,
                       cohort = meta$cohort)

p_post <- ggplot(pca2_df, aes(PC1, PC2, color = platform, shape = cohort)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_classic() +
  labs(title = "Post-ComBat PCA", x = "PC1", y = "PC2")

# Export
ggsave("Figure/Post_ComBat_PCA.svg", p_post, width = 6, height = 5)
cat("✅ Saved: Figure/Post_ComBat_PCA.svg\n")


##### H2  ####

## STEP 1 | Convert Ensembl ID → Gene Symbol


# Install if missing
if (!requireNamespace("biomaRt", quietly = TRUE)) install.packages("biomaRt")
library(biomaRt)
library(dplyr)
library(ggplot2)

# 1) Build Ensembl-to-Symbol mapping
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
map <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), mart = mart)

# 2) Clean Ensembl IDs (remove version suffix, e.g. ENSG00000105352.8 → ENSG00000105352)
rownames_clean <- sub("\\..*", "", rownames(expr_combat))

# 3) Map to gene symbols
gene_symbol <- map$hgnc_symbol[match(rownames_clean, map$ensembl_gene_id)]

# 4) Replace rownames with gene symbols
expr_symbol <- expr_combat[!is.na(gene_symbol) & gene_symbol != "", ]
rownames(expr_symbol) <- gene_symbol[!is.na(gene_symbol) & gene_symbol != ""]

# Optional check
cat("Genes after mapping:", nrow(expr_symbol), "\n")
head(rownames(expr_symbol))

# Replace expr_combat with symbol-level version
expr_combat <- expr_symbol



## STEP 2 | 18-gene CTS composite score

# 1) Define 18-gene signature
cts18_genes <- c(
  "ACER3", "SERPINB1", "HK3", "TDRD9", "NLRC4",
  "PGD", "UBE2H", "METTL9", "STOM", "SNX3",
  "GADD45A", "BTN3A3", "BPGM", "CA1", "SLC4A1",
  "EPB42", "FECH", "GLRX5"
)


# 2) Check which genes are available
cts18_present <- intersect(cts18_genes, rownames(expr_combat))
cat("Number of matched CTS18 genes:", length(cts18_present), "\n")

# 3) Extract expression for those genes
expr_cts18 <- expr_combat[cts18_present, ]

# 4) Z-score per gene (standardize)
expr_cts18_z <- t(scale(t(expr_cts18)))

# 5) Compute per-sample composite score
cts_score <- colMeans(expr_cts18_z, na.rm = TRUE)

# 6) Merge with metadata (make sure meta$sample matches colnames)
meta$CTS_signature_score <- cts_score[meta$sample]

# Ensure cohort factor is set
meta$cohort <- factor(meta$cohort, levels = c("CAP", "Sepsis"))

# 7) Plot violin + boxplot
ggplot(meta, aes(x = cohort, y = CTS_signature_score, fill = cohort)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, outlier.size = 0.5, fill = "white") +
  theme_classic(base_size = 12) +
  labs(title = "18-gene CTS composite score",
       x = NULL, y = "Mean Z-score of 18 genes")

# 8) Wilcoxon test
wilcox.test(CTS_signature_score ~ cohort, data = meta) #p-value = 0.5104 ####



## H5 | PCA of 18 CTS genes #####

library(ggplot2)
library(ggrepel)
library(dplyr)

## Prepare PCA data (already have expr_cts18, studydata_with_cts, dfwb_sel) 

# Run PCA on the 18-gene expression matrix
# expr_cts18: 18 genes × samples (columns are sample names)
pca_res <- prcomp(t(expr_cts18), scale. = TRUE)

# Build sample-level PCA dataframe
pca_df <- data.frame(
  PC1    = pca_res$x[, 1],
  PC2    = pca_res$x[, 2],
  Sample = colnames(expr_cts18)
)

# Load mapping from Sample → ID
dfwb_sel <- read.csv("Original_Data/qns_ID.csv", check.names = FALSE)

# Load CTS info (contains EB_id and CTS)
studydata_with_cts <- read.csv("Original_Data/studydata_with_cts.csv", check.names = FALSE)

# Merge CTS labels to PCA data
# Build sample-level PCA dataframe (CTS 1–3 only already filtered upstream)
# Build sample-level PCA dataframe
pca_df <- data.frame(
  PC1    = pca_res$x[, 1],
  PC2    = pca_res$x[, 2],
  Sample = colnames(expr_cts18)
) %>%
  dplyr::left_join(dplyr::select(dfwb_sel, Sample, ID), by = "Sample") %>%
  dplyr::left_join(dplyr::select(studydata_with_cts, EB_id, CTS),
                   by = c("ID" = "EB_id")) %>%
  dplyr::filter(!is.na(CTS)) %>%   # 去掉NA样本
  dplyr::mutate(
    CTS = factor(CTS, levels = c("1", "2", "3"),
                 labels = c("CTS-1", "CTS-2", "CTS-3"))
  )


#2) Extract gene loadings for arrows
loadings <- as.data.frame(pca_res$rotation[, 1:2])  # PC1 & PC2 loadings
loadings$gene <- rownames(loadings)
colnames(loadings)[1:2] <- c("PC1", "PC2")

# Scale arrow length for visibility
arrow_scale <- 8
loadings <- loadings %>%
  mutate(PC1 = PC1 * arrow_scale,
         PC2 = PC2 * arrow_scale)

## 3) Plot PCA with ellipses + arrows 
p_pca <- ggplot(pca_df, aes(PC1, PC2, color = CTS)) +
  geom_point(size = 2.5, alpha = 0.85) +
  stat_ellipse(type = "norm", level = 0.95, linetype = 2, size = 0.7, alpha = 0.8) +
  # Add gene loading arrows
  geom_segment(
    data = loadings,
    aes(x = 0, y = 0, xend = PC1, yend = PC2),
    arrow = arrow(length = unit(0.2, "cm")),
    color = "gray35", linewidth = 0.5,
    inherit.aes = FALSE
  ) +
  geom_text_repel(
    data = loadings,
    aes(x = PC1, y = PC2, label = gene),
    color = "black", size = 3.2, inherit.aes = FALSE,
    max.overlaps = 100
  ) +
  geom_hline(yintercept = 0, color = "gray85", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "gray85", linewidth = 0.3) +
  scale_color_manual(values = c("CTS-1" = "royalblue",
                                "CTS-2" = "#B2DF8A",
                                "CTS-3" = "orange")) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 11),
    legend.text  = element_text(size = 10),
    plot.title   = element_text(size = 14, face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "PCA of 18 CTS genes (95% ellipse + gene loadings)",
    x = paste0("PC1 (", round(summary(pca_res)$importance[2,1]*100,1), "% variance)"),
    y = paste0("PC2 (", round(summary(pca_res)$importance[2,2]*100,1), "% variance)")
  ) +
  guides(color = guide_legend(title = "CTS subtype",
                              override.aes = list(size = 3)))

# Display on screen
print(p_pca)

## 4) Save figure 
dir.create("Results", showWarnings = FALSE)
ggsave("Results/PCA_18genes_with_ellipse_arrows.svg", p_pca, width = 8, height = 5)
cat("✅ Saved figure: Results/PCA_18genes_with_ellipse_arrows.svg\n")  #####


