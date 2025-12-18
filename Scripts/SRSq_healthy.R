rm(list = ls())

########################################
### 0. Load packages
########################################
library(edgeR)
library(SepstratifieR)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggridges)

########################################
### 1. Load HEALTHY raw counts
########################################
cnt <- read.csv(
  "Original_data/Combined_counts_Amsterdam_my_cohort_allhealthy.csv",
  row.names = 1,
  check.names = FALSE
)

cat("Healthy matrix:", dim(cnt), "\n")

########################################
### 2. Convert raw counts to logCPM
########################################
dge <- DGEList(counts = cnt)
expr_logcpm <- cpm(dge, log = TRUE, prior.count = 1)

# SepstratifieR requires samples x genes
X <- t(expr_logcpm)

########################################
### 3. Run SRS classifier (7 genes)
########################################
pred7 <- stratifyPatients(
  X,
  gene_set = "davenport",
  verbose = TRUE
)

########################################
### 4. Run SRS classifier (19 genes)
########################################
pred19 <- stratifyPatients(
  X,
  gene_set = "extended",
  verbose = TRUE
)

########################################
### 5. Extract results
########################################
res7 <- data.frame(
  Sample = rownames(X),
  SRS_7  = pred7@SRS,
  SRSq_7 = pred7@SRSq,
  stringsAsFactors = FALSE
)

res19 <- data.frame(
  Sample = rownames(X),
  SRS_19  = pred19@SRS,
  SRSq_19 = pred19@SRSq,
  stringsAsFactors = FALSE
)

########################################
### 6. Merge into one dataframe
########################################
out <- res7 %>%
  left_join(res19, by = "Sample")

########################################
### 7. Save results
########################################
write.csv(
  out,
  "Original_data/HEALTHY_SRS_results_7and19.csv",
  row.names = FALSE
)

cat("\n✅ Saved: Original_data/HEALTHY_SRS_results_7and19.csv\n")

########################################
### 8. PCA-aligned sample plots
########################################
p7 <- plotAlignedSamples(pred7) +
  ggtitle("Healthy aligned to SRS reference – 7-gene signature")

p19 <- plotAlignedSamples(pred19) +
  ggtitle("Healthy aligned to SRS reference – 19-gene signature")

p7
p19

########################################
### 9. k-sensitivity analysis
########################################
sens_h <- runSensitivityAnalysis(
  X,
  gene_set = "davenport"
)

print(sens_h)
cat("\n✅ Healthy k-sensitivity complete.\n")

########################################
### 10. FIGURES for healthy
########################################
df <- out
dir.create("Figure", showWarnings = FALSE)

########################################
### 10A. Ridgeline density: SRSq_7 vs SRSq_19
########################################
df_long <- df %>%
  dplyr::select(Sample, SRSq_7, SRSq_19) %>%
  pivot_longer(
    cols = c(SRSq_7, SRSq_19),
    names_to = "Signature",
    values_to = "SRSq"
  )

p_ridge <- ggplot(df_long, aes(x = SRSq, y = Signature, fill = Signature)) +
  geom_density_ridges(alpha = 0.7) +
  labs(
    x = "SRSq",
    y = "Signature",
    title = "Healthy: Distribution of SRSq (7-gene vs 19-gene)"
  ) +
  theme_classic() +
  theme(legend.position = "none")

ggsave(
  "Figure/HEALTHY_SRSq_ridgeline_7vs19.svg",
  p_ridge,
  width = 7,
  height = 5
)

########################################
### 11. Unified SRS color palette (same as PCA)
########################################
srs_cols <- c(
  "SRS1" = "#8B0000",  # dark red
  "SRS2" = "#4F81BD",  # medium blue
  "SRS3" = "#000080"   # dark navy
)

########################################
### 12. PROPORTION PLOT — 7-gene SRS
########################################
prop7 <- df %>%
  group_by(SRS_7) %>%
  summarise(N = n()) %>%
  mutate(Percent = 100 * N / sum(N))

p_prop7 <- ggplot(prop7, aes(x = SRS_7, y = Percent, fill = SRS_7)) +
  geom_col(alpha = 0.9) +
  geom_text(aes(label = sprintf("%.1f%%", Percent)), vjust = -0.3, size = 4) +
  scale_fill_manual(values = srs_cols) +
  labs(
    title = "Healthy: Proportion of SRS groups (7-gene)",
    x = "SRS group",
    y = "Percentage of patients (%)"
  ) +
  theme_classic() +
  theme(legend.position = "none")

ggsave(
  "Figure/HEALTHY_SRS_proportion_7gene.svg",
  p_prop7,
  width = 6,
  height = 5
)

########################################
### 13. PROPORTION PLOT — 19-gene SRS
########################################
prop19 <- df %>%
  group_by(SRS_19) %>%
  summarise(N = n()) %>%
  mutate(Percent = 100 * N / sum(N))

p_prop19 <- ggplot(prop19, aes(x = SRS_19, y = Percent, fill = SRS_19)) +
  geom_col(alpha = 0.9) +
  geom_text(aes(label = sprintf("%.1f%%", Percent)), vjust = -0.3, size = 4) +
  scale_fill_manual(values = srs_cols) +
  labs(
    title = "Healthy: Proportion of SRS groups (19-gene)",
    x = "SRS group",
    y = "Percentage of patients (%)"
  ) +
  theme_classic() +
  theme(legend.position = "none")

ggsave(
  "Figure/HEALTHY_SRS_proportion_19gene.svg",
  p_prop19,
  width = 6,
  height = 5
)

cat("\n✅ All healthy figures saved to folder: Figure/\n")




#cross table
###############################################
### 1. Fix CTS: DO NOT recode 1/2/3
###############################################
cts_df_noninf <- pred %>%
  dplyr::select(Sample, CTS) %>%
  mutate(
    CTS = factor(CTS, levels = c("CTS1","CTS2","CTS3"))
  )

###############################################
### 2. SRS19（无变化）
###############################################
srs_df_noninf <- out %>%
  dplyr::select(Sample, SRS_19) %>% 
  mutate(
    SRS19 = factor(SRS_19,
                   levels = c("SRS1","SRS2","SRS3"))
  )

###############################################
### 3. Merge CTS + SRS19
###############################################
cts_srs_noninf <- inner_join(
  cts_df_noninf,
  srs_df_noninf,
  by = "Sample"
)

###############################################
### 4. Cross-table
###############################################
cross_tab_noninf <- table(cts_srs_noninf$CTS, cts_srs_noninf$SRS19)
cross_tab_noninf

prop_noninf <- round(prop.table(cross_tab_noninf, 1) * 100, 1)
prop_noninf

###############################################
### 5. Plot
###############################################
library(ggplot2)

srs_cols <- c("SRS1"="#8B0000", 
              "SRS2"="#4F81BD", 
              "SRS3"="#000080")

p_noninf <- ggplot(cts_srs_noninf, aes(x = CTS, fill = SRS19)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = srs_cols) +
  labs(
    title = "Non-infection: Proportion of SRS19 within each CTS",
    x = "CTS subtype",
    y = "Percent"
  ) +
  theme_classic(base_size = 14)

p_noninf
print(p_noninf)
dev.off()   # 如果图形设备卡住，运行这个
p_noninf





############################################################
### 1. CTS dataframe
############################################################
cts_df <- pred %>%
  dplyr::select(Sample, CTS) %>%
  mutate(
    CTS = factor(CTS,
                 levels = c("CTS1","CTS2","CTS3"))
  )

############################################################
### 2. SRS7 dataframe
############################################################
srs7_df <- out %>%
  dplyr::select(Sample, SRS_7) %>%
  mutate(
    SRS7 = factor(SRS_7,
                  levels = c("SRS1","SRS2","SRS3"))
  )

############################################################
### 3. Merge CTS + SRS7
############################################################
cts_srs7 <- inner_join(cts_df, srs7_df, by = "Sample")

head(cts_srs7)
table(cts_srs7$CTS)
table(cts_srs7$SRS7)

############################################################
### 4. Stacked proportion plot
############################################################
library(ggplot2)

srs_cols <- c(
  "SRS1" = "#8B0000",
  "SRS2" = "#4F81BD",
  "SRS3" = "#000080"
)

p_srs7 <- ggplot(cts_srs7, aes(x = CTS, fill = SRS7)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = srs_cols) +
  labs(
    title = "Non-infection: Proportion of SRS7 within each CTS",
    x = "CTS subtype",
    y = "Percent"
  ) +
  theme_classic(base_size = 14)

p_srs7

p_srs7
print(p_srs7)
dev.off()   # 如果图形设备卡住，运行这个
p_srs7
