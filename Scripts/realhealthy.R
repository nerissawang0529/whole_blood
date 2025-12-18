rm(list = ls())

## Healthy volunteers Nerissa U219 / HTA data 

# -----------------------------------------------------------------------------#
# ------------ Compare healthies for coconut ----------------------------------#
# -----------------------------------------------------------------------------#

# Load packages 
library(dplyr)
library(tableone)
library(lubridate)
library(DESeq2)
library(rio)
library(GEOquery)
library(ggplot2)



#CTS

# U219 ------------------------------------------------------------------------
load("Original_data/MARS_U219_correct_age.eset")

# inspection of healthy volunteers=> 42 available 
U219_hv <- pData(U219_correct_ages) %>%
  dplyr::filter(health == "healthy") %>%
  dplyr::select(MARSID, Patient_gender, age)

# filter out healthies from entire eset 
healthy_eset_U219 <- U219_correct_ages[, U219_correct_ages$health == "healthy"]
dim(healthy_eset_U219) # correct, 42 samples 

expr_mat <- exprs(healthy_eset_U219)

library(hgu219.db)

# Get mapping
ens_map <- AnnotationDbi::select(
  hgu219.db,
  keys = rownames(expr_mat),
  columns = c("ENSEMBL"),
  keytype = "PROBEID"
)

# Remove probes without ENSEMBL mapping
ens_map <- ens_map[!is.na(ens_map$ENSEMBL), ]

# Collapse to unique ENSEMBL (if multiple probes map to same ENSEMBL, keep mean)
library(dplyr)

expr_ens <- expr_mat %>%
  as.data.frame() %>%
  mutate(PROBEID = rownames(expr_mat)) %>%
  inner_join(ens_map, by = "PROBEID") %>%
  group_by(ENSEMBL) %>%
  summarise(across(where(is.numeric), mean)) %>%
  tibble::column_to_rownames("ENSEMBL")

library(ConsensusTranscriptomicSubtype)

result <- run_subtype_classifier(new_expr_data = expr_ens)

result$predictions
result$probabilities


library(ComplexHeatmap)
library(circlize)
library(org.Hs.eg.db)
library(AnnotationDbi)

## 1. Extract expression + CTS
mat  <- result$expression_corrected$new
pred <- result$predictions

pred$CTS <- factor(pred$CTS,
                   levels = c("1","2","3"),
                   labels = c("CTS1","CTS2","CTS3"))

## 2. ENSG → SYMBOL
ens_ids <- sub("\\..*$", "", rownames(mat))

gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys     = ens_ids,
  keytype  = "ENSEMBL",
  column   = "SYMBOL",
  multiVals = "first"
)

rownames(mat) <- gene_symbols
mat <- mat[!is.na(rownames(mat)) & rownames(mat) != "", , drop = FALSE]

## 3. Keep genes in fixed order
gene_order <- c(
  "ACER3","SERPINB1","HK3","TDRD9","NLRC4","PGD","UBE2H",
  "METTL9","STOM","SNX3","GADD45A","BTN3A3","BPGM",
  "CA1","SLC4A1","EPB42","FECH","GLRX5"
)

mat2 <- mat[gene_order, , drop = FALSE]

## 4. Row Z-score
mat_scaled <- t(scale(t(mat2)))

## 🔹 5. Reorder samples by CTS: CTS1 → CTS2 → CTS3
ord <- order(pred$CTS)          # factor order handles 1,2,3
mat_scaled <- mat_scaled[, ord]
pred       <- pred[ord, ]

## 6. Annotation + colors
cts_cols <- c("CTS1"="royalblue",
              "CTS2"="#B2DF8A",
              "CTS3"="orange")

top_annot <- HeatmapAnnotation(
  CTS = pred$CTS,
  col  = list(CTS = cts_cols),
  show_annotation_name = FALSE
)

col_fun <- colorRamp2(
  c(-3, 0, 3),
  c("#4B4B9A","white","#C47E32")
)

## 7. Draw heatmap in RStudio
Heatmap(
  mat_scaled,
  name = "Z-score",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 9),
  top_annotation = top_annot
)

## Count per CTS
cts_counts <- table(pred$CTS)
cts_title  <- paste0(
  "CTS1 (n=", cts_counts["CTS1"], ")    ",
  "CTS2 (n=", cts_counts["CTS2"], ")    ",
  "CTS3 (n=", cts_counts["CTS3"], ")"
)

## ============================
## 7. Draw heatmap in RStudio with CTS counts
## ============================

# Count CTS1/2/3
cts_counts <- table(pred$CTS)
cts_title  <- paste0(
  "CTS1 (n=", cts_counts["CTS1"], ")    ",
  "CTS2 (n=", cts_counts["CTS2"], ")    ",
  "CTS3 (n=", cts_counts["CTS3"], ")"
)

Heatmap(
  mat_scaled,
  name = "Z-score",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 9),
  top_annotation = top_annot,
  column_title = cts_title,
  column_title_gp = gpar(fontsize = 12, fontface = "bold")
)




#select the patients <= 50 years old
# Select HV < 50 years
hv_under50 <- U219_hv %>% 
  filter(age < 50)

hv_under50_ids <- hv_under50$MARSID
length(hv_under50_ids)

library(stringr)

# Extract hvID from column names (everything after last "_")
expr_ids <- str_extract(colnames(expr_ens), "hv\\d+")

# Keep only columns matching hv_under50_ids
expr_ens_under50 <- expr_ens[, expr_ids %in% hv_under50_ids]
dim(expr_ens_under50)


result_under50 <- run_subtype_classifier(new_expr_data = expr_ens_under50)

result_under50$predictions

library(ComplexHeatmap)
library(circlize)
library(org.Hs.eg.db)
library(AnnotationDbi)

## 1. Extract expression + CTS
mat  <- result_under50$expression_corrected$new
pred <- result_under50$predictions

pred$CTS <- factor(pred$CTS,
                   levels = c("1","2","3"),
                   labels = c("CTS1","CTS2","CTS3"))

## 2. ENSG → SYMBOL
ens_ids <- sub("\\..*$", "", rownames(mat))

gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys     = ens_ids,
  keytype  = "ENSEMBL",
  column   = "SYMBOL",
  multiVals = "first"
)

rownames(mat) <- gene_symbols
mat <- mat[!is.na(rownames(mat)) & rownames(mat) != "", , drop = FALSE]

## 3. Keep classifier genes
gene_order <- c(
  "ACER3","SERPINB1","HK3","TDRD9","NLRC4","PGD","UBE2H",
  "METTL9","STOM","SNX3","GADD45A","BTN3A3","BPGM",
  "CA1","SLC4A1","EPB42","FECH","GLRX5"
)
mat2 <- mat[gene_order, , drop = FALSE]

## 4. Z-score
mat_scaled <- t(scale(t(mat2)))

## 5. Order samples CTS1 → CTS2 → CTS3
ord <- order(pred$CTS)
mat_scaled <- mat_scaled[, ord]
pred       <- pred[ord, ]

## 6. CTS colors
cts_cols <- c("CTS1"="royalblue",
              "CTS2"="#B2DF8A",
              "CTS3"="orange")

top_annot <- HeatmapAnnotation(
  CTS = pred$CTS,
  col  = list(CTS = cts_cols),
  show_annotation_name = FALSE
)

col_fun <- colorRamp2(c(-3,0,3),
                      c("#4B4B9A","white","#C47E32"))

## 7. Add CTS counts
cts_counts <- table(pred$CTS)
cts_title  <- paste0(
  "CTS1 (n=", cts_counts["CTS1"], ")    ",
  "CTS2 (n=", cts_counts["CTS2"], ")    ",
  "CTS3 (n=", cts_counts["CTS3"], ")"
)

## 8. Draw heatmap
Heatmap(
  mat_scaled,
  name = "Z-score",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 9),
  top_annotation = top_annot,
  column_title = cts_title,
  column_title_gp = gpar(fontsize = 12, fontface = "bold")
)



#SRS####

########################################
### Load packages
########################################
library(SepstratifieR)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggridges)
library(stringr)

########################################
### Input variables from your workflow:
### expr_ens = ENSEMBL × 42 HV microarray expression
### U219_hv = demographics (MARSID, gender, age)
########################################

expr_all <- expr_ens   # rename for clarity
cat("U219 HV matrix:", dim(expr_all), "\n")

# Convert to samples × genes (SepstratifieR format)
X_all <- t(expr_all)

########################################
### OFFICIAL SRS COLORS
########################################
srs_cols <- c(
  "SRS1" = "#8B0000",   # dark red
  "SRS2" = "#4F81BD",   # bright blue
  "SRS3" = "#000080"    # dark navy
)

########################################
### 1️⃣ RUN SRS — ALL 42 HV
########################################

pred7_all  <- stratifyPatients(X_all, gene_set = "davenport")
pred19_all <- stratifyPatients(X_all, gene_set = "extended")

# Combine results
res_all <- data.frame(
  Sample = rownames(X_all),
  SRS7   = pred7_all@SRS,
  SRSq7  = pred7_all@SRSq,
  SRS19  = pred19_all@SRS,
  SRSq19 = pred19_all@SRSq
)

dir.create("Original_data", showWarnings = FALSE)

write.csv(res_all,
          "Original_data/U219_HEALTHY_SRS_results_7and19.csv",
          row.names = FALSE)

cat("\n✅ Saved: U219_HEALTHY_SRS_results_7and19.csv\n")

########################################
### 2️⃣ PCA ALIGNMENT — ALL HV (CORRECT COLORS)
########################################

p7_all <- plotAlignedSamples(pred7_all) +
  scale_color_manual(values = srs_cols, name = "SRS") +
  ggtitle("U219 Healthy – SRS 7-gene alignment")

p19_all <- plotAlignedSamples(pred19_all) +
  scale_color_manual(values = srs_cols, name = "SRS") +
  ggtitle("U219 Healthy – SRS 19-gene alignment")

p7_all
p19_all

########################################
### 3️⃣ PROPORTION PLOTS — ALL HV
########################################

dir.create("Figure", showWarnings = FALSE)

# 7-gene
prop7 <- res_all %>%
  group_by(SRS7) %>%
  summarise(N = n()) %>%
  mutate(Percent = 100*N/sum(N))

p_prop7 <- ggplot(prop7, aes(x = SRS7, y = Percent, fill = SRS7)) +
  geom_col() +
  geom_text(aes(label = sprintf("%.1f%%", Percent)), vjust = -0.3) +
  scale_fill_manual(values = srs_cols) +
  labs(title="U219 Healthy – SRS Proportion (7-gene)") +
  theme_classic()

ggsave("Figure/U219_SRS_proportion_7gene.svg",
       p_prop7, width = 6, height = 5)

# 19-gene
prop19 <- res_all %>%
  group_by(SRS19) %>%
  summarise(N = n()) %>%
  mutate(Percent = 100*N/sum(N))

p_prop19 <- ggplot(prop19, aes(x = SRS19, y = Percent, fill = SRS19)) +
  geom_col() +
  geom_text(aes(label = sprintf("%.1f%%", Percent)), vjust = -0.3) +
  scale_fill_manual(values = srs_cols) +
  labs(title="U219 Healthy – SRS Proportion (19-gene)") +
  theme_classic()

ggsave("Figure/U219_SRS_proportion_19gene.svg",
       p_prop19, width = 6, height = 5)

########################################
### 4️⃣ RUN SRS — AGE ≤ 50 ONLY
########################################

hv_under50 <- U219_hv %>% filter(age <= 50)
hv_under50_ids <- hv_under50$MARSID

expr_ids <- str_extract(colnames(expr_ens), "hv\\d+")
expr_u50 <- expr_ens[, expr_ids %in% hv_under50_ids]

X_u50 <- t(expr_u50)

pred7_u50  <- stratifyPatients(X_u50, gene_set = "davenport")
pred19_u50 <- stratifyPatients(X_u50, gene_set = "extended")

res_u50 <- data.frame(
  Sample = rownames(X_u50),
  SRS7   = pred7_u50@SRS,
  SRSq7  = pred7_u50@SRSq,
  SRS19  = pred19_u50@SRS,
  SRSq19 = pred19_u50@SRSq
)

write.csv(res_u50,
          "Original_data/U219_HEALTHY_ageUnder50_SRS_results.csv",
          row.names = FALSE)

cat("\n✅ Saved: U219_HEALTHY_ageUnder50_SRS_results.csv\n")

########################################
### 5️⃣ PCA ALIGNMENT — AGE ≤ 50 (CORRECT COLORS)
########################################

########################################
### 5️⃣ PCA ALIGNMENT — AGE ≤ 50 (CORRECT COLORS)
########################################

p7_u50 <- plotAlignedSamples(pred7_u50) +
  scale_color_manual(values = srs_cols, name = "SRS") +
  ggtitle("U219 Healthy (≤50) – SRS 7-gene alignment")

# Flip PC1 horizontally so orientation matches the other plots
p19_u50 <- plotAlignedSamples(pred19_u50) +
  scale_color_manual(values = srs_cols, name = "SRS") +
  ggtitle("U219 Healthy (≤50) – SRS 19-gene alignment") +
  scale_x_reverse()

p7_u50
p19_u50


########################################
cat("\n\n🎉 ALL DONE!\n")
cat("• SRS for all HV saved\n")
cat("• SRS for ≤50 saved\n")
cat("• PCA figures displayed (correct colors)\n")
cat("• Proportion figures saved in Figure/\n")


########################################
### 6️⃣ PROPORTION PLOTS — AGE ≤ 50 HV
########################################

# 7-gene proportions (≤50)
prop7_u50 <- res_u50 %>%
  group_by(SRS7) %>%
  summarise(N = n()) %>%
  mutate(Percent = 100*N/sum(N))

p_prop7_u50 <- ggplot(prop7_u50, aes(x = SRS7, y = Percent, fill = SRS7)) +
  geom_col() +
  geom_text(aes(label = sprintf("%.1f%%", Percent)), vjust = -0.3) +
  scale_fill_manual(values = srs_cols) +
  labs(title = "U219 Healthy (≤50) – SRS Proportion (7-gene)") +
  theme_classic()

ggsave("Figure/U219_SRS_proportion_7gene_under50.svg",
       p_prop7_u50, width = 6, height = 5)


# 19-gene proportions (≤50)
prop19_u50 <- res_u50 %>%
  group_by(SRS19) %>%
  summarise(N = n()) %>%
  mutate(Percent = 100*N/sum(N))

p_prop19_u50 <- ggplot(prop19_u50, aes(x = SRS19, y = Percent, fill = SRS19)) +
  geom_col() +
  geom_text(aes(label = sprintf("%.1f%%", Percent)), vjust = -0.3) +
  scale_fill_manual(values = srs_cols) +
  labs(title = "U219 Healthy (≤50) – SRS Proportion (19-gene)") +
  theme_classic()

ggsave("Figure/U219_SRS_proportion_19gene_under50.svg",
       p_prop19_u50, width = 6, height = 5)

cat("\n📊 Saved: proportion plots for ≤50\n")


########################################
### 6️⃣ PROPORTION PLOTS — AGE ≤ 50 HV
########################################

# 7-gene proportions (≤50)
prop7_u50 <- res_u50 %>%
  group_by(SRS7) %>%
  summarise(N = n()) %>%
  mutate(Percent = 100*N/sum(N))

p_prop7_u50 <- ggplot(prop7_u50, aes(x = SRS7, y = Percent, fill = SRS7)) +
  geom_col() +
  geom_text(aes(label = sprintf("%.1f%%", Percent)), vjust = -0.3) +
  scale_fill_manual(values = srs_cols) +
  labs(title = "U219 Healthy (≤50) – SRS Proportion (7-gene)") +
  theme_classic()

ggsave("Figure/U219_SRS_proportion_7gene_under50.svg",
       p_prop7_u50, width = 6, height = 5)


# 19-gene proportions (≤50)
prop19_u50 <- res_u50 %>%
  group_by(SRS19) %>%
  summarise(N = n()) %>%
  mutate(Percent = 100*N/sum(N))

p_prop19_u50 <- ggplot(prop19_u50, aes(x = SRS19, y = Percent, fill = SRS19)) +
  geom_col() +
  geom_text(aes(label = sprintf("%.1f%%", Percent)), vjust = -0.3) +
  scale_fill_manual(values = srs_cols) +
  labs(title = "U219 Healthy (≤50) – SRS Proportion (19-gene)") +
  theme_classic()

ggsave("Figure/U219_SRS_proportion_19gene_under50.svg",
       p_prop19_u50, width = 6, height = 5)

cat("\n📊 Saved: proportion plots for ≤50\n")








########################################
### CHECK SRS19 GENE COVERAGE IN U219
########################################

library(org.Hs.eg.db)
library(AnnotationDbi)

## 1. SRS19 extended signature (Ensembl IDs)
srs19_ensembl <- c(
  "ENSG00000144659", # SLC25A38
  "ENSG00000103423", # DNAJA3
  "ENSG00000135372", # NAT10
  "ENSG00000079134", # THOC1
  "ENSG00000135972", # MRPS9
  "ENSG00000087157", # PGS1
  "ENSG00000165006", # UBAP1
  "ENSG00000111667", # USP5
  "ENSG00000182670", # TTC3
  "ENSG00000097033", # SH3GLB1
  "ENSG00000165733", # BMS1
  "ENSG00000103264", # FBXO31
  "ENSG00000152219", # ARL14EP
  "ENSG00000100814", # CCNB1IP1
  "ENSG00000127334", # DYRK2
  "ENSG00000131355", # ADGRE3
  "ENSG00000137337", # MDC1
  "ENSG00000156414", # TDRD9
  "ENSG00000115085"  # ZAP70
)

cat("SRS19 genes (Ensembl):", length(srs19_ensembl), "\n")

## 2. Make sure expr_ens rownames are Ensembl IDs without version
u219_ens <- sub("\\..*$", "", rownames(expr_ens))

## 3. Present vs missing (by Ensembl ID)
present_ids <- srs19_ensembl[srs19_ensembl %in% u219_ens]
missing_ids <- srs19_ensembl[!srs19_ensembl %in% u219_ens]

## 4. Map to symbols for readability
present_symbols <- mapIds(
  org.Hs.eg.db,
  keys     = present_ids,
  keytype  = "ENSEMBL",
  column   = "SYMBOL",
  multiVals = "first"
)

missing_symbols <- mapIds(
  org.Hs.eg.db,
  keys     = missing_ids,
  keytype  = "ENSEMBL",
  column   = "SYMBOL",
  multiVals = "first"
)

## 5. Print nicely
cat("\n================ PRESENT SRS19 GENES IN U219 ================\n")
print(data.frame(
  ENSEMBL = present_ids,
  SYMBOL  = present_symbols
))

cat("\n================ MISSING SRS19 GENES IN U219 ================\n")
print(data.frame(
  ENSEMBL = missing_ids,
  SYMBOL  = missing_symbols
))
cat("=============================================================\n")



########################################
### Combine CTS + SRS table
########################################

cts_df <- data.frame(
  Sample = colnames(mat_scaled),   # from CTS heatmap
  CTS    = pred$CTS
)

srs_df <- res_all %>% 
  dplyr::select(Sample, SRS19, SRSq19)

# Merge
cts_srs <- inner_join(cts_df, srs_df, by = "Sample")

head(cts_srs)

########################################
### CTS × SRS cross table
########################################

cross_tab <- table(cts_srs$CTS, cts_srs$SRS19)
cross_tab

prop_table <- prop.table(cross_tab, margin = 1) * 100
round(prop_table, 1)

########################################
### Barplot of SRS in each CTS
########################################

library(ggplot2)

ggplot(cts_srs, aes(x = CTS, fill = SRS19)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = srs_cols) +
  labs(
    title = "Proportion of SRS19 classes in each CTS",
    y = "Percent"
  ) +
  theme_classic()

########################################
### SRS7 × CTS (同 SRS19 风格)
########################################

# 假设你已有：
# cts_df = data.frame(Sample, CTS)     # 来自 CTS heatmap
# res_all = data.frame(Sample, SRS7...) # SRS7/SRS19 结果

srs7_df <- res_all %>% 
  dplyr::select(Sample, SRS7)

cts_srs7 <- inner_join(cts_df, srs7_df, by = "Sample")

# 只画 SRS1/SRS2/SRS3
cts_srs7$SRS7 <- factor(cts_srs7$SRS7, levels=c("SRS1","SRS2","SRS3"))

# 与你之前完全相同的配色
srs_cols <- c(
  "SRS1"="#8B0000",
  "SRS2"="#4F81BD",
  "SRS3"="#000080"
)

p_srs7 <- ggplot(cts_srs7, aes(x = CTS, fill = SRS7)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = srs_cols) +
  labs(
    title = "Proportion of SRS7 classes in each CTS",
    y = "Percent"
  ) +
  theme_classic()

p_srs7
###############################################
### 1. Get ≤50 CTS results
###############################################

# CTS (≤50) obtained earlier
mat_u50  <- result_under50$expression_corrected$new
pred_u50 <- result_under50$predictions

pred_u50$CTS <- factor(pred_u50$CTS,
                       levels = c("1","2","3"),
                       labels = c("CTS1","CTS2","CTS3"))

cts_df_u50 <- data.frame(
  Sample = pred_u50$Sample,
  CTS    = pred_u50$CTS
)

###############################################
### 2. Get ≤50 SRS results
###############################################

# Already computed:
# res_u50$Sample, res_u50$SRS7, res_u50$SRS19

srs_df_u50 <- res_u50 %>%
  dplyr::select(Sample, SRS7, SRS19)

###############################################
### 3. Merge CTS + SRS for ≤50
###############################################

cts_srs_u50 <- inner_join(cts_df_u50, srs_df_u50, by = "Sample")

###############################################
### 4. Official SRS colors
###############################################

srs_cols <- c(
  "SRS1" = "#8B0000",
  "SRS2" = "#4F81BD",
  "SRS3" = "#000080"
)

###############################################
### 5. Plot — SRS7 × CTS (≤50)
###############################################

library(ggplot2)

p_srs7_u50 <- ggplot(cts_srs_u50, aes(x = CTS, fill = SRS7)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = srs_cols) +
  labs(
    title = "Healthy ≤50: Proportion of SRS7 within each CTS",
    x = "CTS subtype",
    y = "Percent"
  ) +
  theme_classic(base_size = 14)

p_srs7_u50


###############################################
### 6. Plot — SRS19 × CTS (≤50)
###############################################

p_srs19_u50 <- ggplot(cts_srs_u50, aes(x = CTS, fill = SRS19)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = srs_cols) +
  labs(
    title = "Healthy ≤50: Proportion of SRS19 within each CTS",
    x = "CTS subtype",
    y = "Percent"
  ) +
  theme_classic(base_size = 14)

p_srs19_u50
