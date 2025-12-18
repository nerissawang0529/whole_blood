# this is for CAP
rm(list = ls())

#######################################
### 0. Load packages
#######################################
library(edgeR)
library(SepstratifieR)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggridges)

#######################################
### 1. Load CAP raw counts (genes × samples)
#######################################
cnt <- read.csv("Original_data/Combined_counts_Amsterdam_my_cohort_allpatients.csv",
                row.names = 1, check.names = FALSE)

dim(cnt)
head(cnt[, 1:5])

#######################################
### 2. Convert to logCPM
#######################################
dge <- DGEList(counts = cnt)
expr_logcpm <- cpm(dge, log = TRUE, prior.count = 1)

# SepstratifieR expects samples × genes
X <- t(expr_logcpm)

#######################################
### 3. Run SRS using 7-gene signature
#######################################
pred7 <- stratifyPatients(
  X,
  gene_set = "davenport",   # 7-gene Davenport signature
  verbose = TRUE
)

#######################################
### 4. Run SRS using 19-gene signature
#######################################
pred19 <- stratifyPatients(
  X,
  gene_set = "extended",    # 19-gene Extended signature
  verbose = TRUE
)

#######################################
### 5. Extract SRS + SRSq
#######################################
res7 <- data.frame(
  Sample = rownames(X),
  SRS_7  = pred7@SRS,       # SRS1 / SRS2 / SRS3 (7-gene)
  SRSq_7 = pred7@SRSq,      # continuous SRSq (7-gene)
  stringsAsFactors = FALSE
)

res19 <- data.frame(
  Sample = rownames(X),
  SRS_19  = pred19@SRS,     # SRS1 / SRS2 / SRS3 (19-gene)
  SRSq_19 = pred19@SRSq,    # continuous SRSq (19-gene)
  stringsAsFactors = FALSE
)

#######################################
### 6. Merge results
#######################################
out <- res7 %>%
  dplyr::left_join(res19, by = "Sample")

#######################################
### 7. Save to Original_data/
#######################################
dir.create("Original_data", showWarnings = FALSE)

write.csv(
  out,
  "Original_data/CAP_SRS_results_7and19.csv",
  row.names = FALSE
)

cat("\n✅ Results saved to: Original_data/CAP_SRS_results_7and19.csv\n")

#######################################
### 8. QC: Plot alignment to reference (PCA-like)
#######################################
p7 <- plotAlignedSamples(pred7) +
  ggtitle("CAP aligned to SRS reference – 7-gene signature")

p19 <- plotAlignedSamples(pred19) +
  ggtitle("CAP aligned to SRS reference – 19-gene signature")

p7
p19

#######################################
### 9. k-sensitivity analysis (run ONCE)
#######################################
sens <- runSensitivityAnalysis(
  X,
  gene_set = "davenport"    # use the 7-gene set
)

print(sens)
cat("\n✅ k-sensitivity (robustness of mNN alignment) completed.\n")

############################################################
### 10. Read SRS results for plots
############################################################
df <- read.csv("Original_data/CAP_SRS_results_7and19.csv")

# Create figure folder
dir.create("Figure", showWarnings = FALSE)

############################################################
### 11. Ridgeline plot of SRSq distributions
############################################################
df_long <- df %>%
  dplyr::select(Sample, SRSq_7, SRSq_19) %>%
  tidyr::pivot_longer(
    cols      = c(SRSq_7, SRSq_19),
    names_to  = "Signature",
    values_to = "SRSq"
  )

p_ridge <- ggplot(df_long, aes(x = SRSq, y = Signature, fill = Signature)) +
  geom_density_ridges(alpha = 0.7) +
  labs(
    x = "SRSq",
    y = "Signature",
    title = "Distribution of SRSq for 7-gene and 19-gene signatures"
  ) +
  theme_classic() +
  theme(legend.position = "none")

ggsave(
  filename = "Figure/SRSq_ridgeline_7vs19.svg",
  plot     = p_ridge,
  width    = 7,
  height   = 5
)

############################################################
### 12. Unified SRS color palette (match PCA)
############################################################
srs_cols <- c(
  "SRS1" = "#8B0000",  # dark red
  "SRS2" = "#4F81BD",  # medium blue
  "SRS3" = "#000080"   # dark navy
)

############################################################
### 13. Proportion plot for 7-gene SRS (colors = SRS palette)
############################################################
prop7 <- df %>%
  dplyr::group_by(SRS_7) %>%
  dplyr::summarise(N = dplyr::n()) %>%
  dplyr::mutate(Percent = 100 * N / sum(N))

p_prop_7 <- ggplot(prop7, aes(x = SRS_7, y = Percent, fill = SRS_7)) +
  geom_col(alpha = 0.9) +
  geom_text(aes(label = sprintf("%.1f%%", Percent)),
            vjust = -0.3, size = 4) +
  scale_fill_manual(values = srs_cols) +
  labs(
    title = "Proportion of patients in each SRS group (7-gene signature)",
    x     = "SRS group",
    y     = "Percentage of patients (%)"
  ) +
  theme_classic() +
  theme(legend.position = "none")

ggsave(
  filename = "Figure/SRS_proportion_7gene.svg",
  plot     = p_prop_7,
  width    = 6,
  height   = 5
)

############################################################
### 14. Proportion plot for 19-gene SRS (colors = SRS palette)
############################################################
prop19 <- df %>%
  dplyr::group_by(SRS_19) %>%
  dplyr::summarise(N = dplyr::n()) %>%
  dplyr::mutate(Percent = 100 * N / sum(N))

p_prop_19 <- ggplot(prop19, aes(x = SRS_19, y = Percent, fill = SRS_19)) +
  geom_col(alpha = 0.9) +
  geom_text(aes(label = sprintf("%.1f%%", Percent)),
            vjust = -0.3, size = 4) +
  scale_fill_manual(values = srs_cols) +
  labs(
    title = "Proportion of patients in each SRS group (19-gene signature)",
    x     = "SRS group",
    y     = "Percentage of patients (%)"
  ) +
  theme_classic() +
  theme(legend.position = "none")

ggsave(
  filename = "Figure/SRS_proportion_19gene.svg",
  plot     = p_prop_19,
  width    = 6,
  height   = 5
)

cat("\n✅ Ridgeline + SRS proportion plots exported to: Figure/\n")
