rm(list = ls ())

# ==========================
# Ready for RNA-seq data — HV cohort
# ==========================

#make CTS for healthy#######
## ===== 1. Read clinical data and identify HV patients =====
eb <- read.csv("Original_data/studydata_all_EB.csv")

# Identify Healthy Volunteers (HV)
hv_patients <- eb[grepl("HV", eb$group, ignore.case = TRUE), ]

# Extract unique Participant IDs for HV
hv_ids <- unique(hv_patients$EB_id)

length(hv_ids)
# ✅ These are all HV patients from the ELDER-BIOME cohort


# ==========================
# Ready for RNA-seq data — HV cohort (BO NL... SG#### pattern)
# ==========================

## 1) Read key and counts
key       <- read.csv("Original_data/SK_Amsterdam.csv")
count_mat <- read.csv("Original_data/Combined_counts_Amsterdam.csv")

# Gene IDs as rownames (strip version suffix after ".")
row.names(count_mat) <- sub("\\..*$", "", count_mat$X)
count_mat$X <- NULL

## 2) Identify HV samples by Sample.ID pattern
# Pattern: "BO NL<digits>.<digits>.<digits> SG<digits>"
hv_pat <- "^BO\\sNL\\d+\\.\\d+\\.\\d+\\sSG\\d+$"
Elder_hv <- subset(key, grepl(hv_pat, Sample.ID) & RNAseq == 1)

# (Optional) quick check
cat("HV samples detected:", nrow(Elder_hv), "\n")
if (!all(Elder_hv$Sample %in% colnames(count_mat))) {
  missing_cols <- setdiff(Elder_hv$Sample, colnames(count_mat))
  cat("⚠️ Missing in count matrix:", length(missing_cols), "columns\n")
}

## 3) Extract RNA-seq expression data for these HV samples
count_hv <- count_mat[, colnames(count_mat) %in% Elder_hv$Sample, drop = FALSE]
cat("HV expression matrix dimensions:", dim(count_hv)[1], "genes x", dim(count_hv)[2], "samples\n")

## 4) Exports
write.csv(Elder_hv, "Original_data/qns_ID_HV.csv", row.names = FALSE)
write.csv(count_hv, "Original_data/Combined_counts_Amsterdam_my_cohort_allhealthy.csv", row.names = TRUE)

cat("✅ Exported:\n",
    "- Original_data/qns_ID_HV.csv\n",
    "- Original_data/Combined_counts_Amsterdam_my_cohort_allhealthy.csv\n")




## ==========================
## 1) Load data
## ==========================
library(tidyverse)
library(edgeR)
library(ConsensusTranscriptomicSubtype)

# Expression counts for healthy samples
cnt <- read.csv("Original_data/Combined_counts_Amsterdam_my_cohort_allhealthy.csv", check.names = FALSE)

# Set gene IDs as rownames (remove Ensembl version suffix)
ens_raw   <- as.character(cnt[[1]])
ens_clean <- sub("\\..*$", "", ens_raw)
cnt[[1]]  <- NULL
rownames(cnt) <- ens_clean

# Convert to numeric matrix
cnt_mat <- as.matrix(cnt)
storage.mode(cnt_mat) <- "numeric"

## --- logCPM transform
dge <- DGEList(counts = cnt_mat)
expr_logcpm <- cpm(dge, log = TRUE, prior.count = 1)

## ==========================
## 2) Prepare for CTS
## ==========================
# Load CTS reference objects
load("Original_data/exp_core_g.rda")
load("Original_data/core_samples.rda")

# Check missing genes
missing <- setdiff(rownames(exp_core_g), rownames(expr_logcpm))
if (length(missing)) stop(paste("Missing CTS genes:", paste(missing, collapse=", ")))

# Align rows
new_expr <- expr_logcpm[rownames(exp_core_g), , drop = FALSE]

## ==========================
## 3) Run CTS classification
## ==========================
dir.create("Figure", showWarnings = FALSE)

res_healthy <- run_subtype_classifier(
  new_expr_data   = new_expr,
  exp_core_g      = exp_core_g,
  core_samples    = core_samples,
  heatmap_file    = "Figure/cts_yourcohort_healthy_heatmap.pdf",
  silhouette_file = "Figure/cts_yourcohort_healthy_silhouette.pdf"
)

# Save predictions
write.csv(res_healthy$predictions, "Figure/cts_yourcohort_healthy_predictions.csv", row.names = FALSE)

## ==========================
## 4) Plot CTS distribution
## ==========================
pred <- res_healthy$predictions
pred$CTS <- factor(pred$CTS, levels = c("1","2","3"))

cts_counts <- as.data.frame(table(pred$CTS))
names(cts_counts) <- c("CTS", "N")

p <- ggplot(cts_counts, aes(x = CTS, y = N)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = N), vjust = -0.3, size = 4) +
  labs(
    title = "CTS subtype counts (Healthy cohort)",
    x = "CTS subtype", y = "Number of samples"
  ) +
  theme_classic(base_size = 12)

ggsave("Figure/cts_yourcohort_subtype_healthy_ggplot.pdf", p, width = 5, height = 4)
cat("✅ CTS classification for healthy cohort completed.\n",
    "Results saved in 'Figure/' folder.\n")

# ======================================
# BEAUTIFIED CTS HEATMAP (Nature style)
# ======================================
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)

# 1️⃣ Load predictions and matrix
pred <- read.csv("Figure/cts_yourcohort_healthy_predictions.csv", check.names = FALSE)
pred$CTS <- factor(pred$CTS, levels = c("1","2","3"), labels = c("CTS1","CTS2","CTS3"))

mat <- res_healthy$expression_corrected$new[, pred$Sample, drop = FALSE]

# 2️⃣ Convert Ensembl → Symbol
if (any(grepl("^ENSG", rownames(mat)))) {
  ens_ids <- sub("\\..*$", "", rownames(mat))
  gene_symbols <- mapIds(org.Hs.eg.db,
                         keys = ens_ids,
                         column = "SYMBOL",
                         keytype = "ENSEMBL",
                         multiVals = "first")
  rownames(mat) <- gene_symbols
}

mat <- mat[!is.na(rownames(mat)) & rownames(mat) != "", , drop = FALSE]

# 3️⃣ Keep only the 18 classifier genes, in fixed order
gene_order <- c("ACER3","SERPINB1","HK3","TDRD9","NLRC4","PGD","UBE2H",
                "METTL9","STOM","SNX3","GADD45A","BTN3A3","BPGM",
                "CA1","SLC4A1","EPB42","FECH","GLRX5")
mat <- mat[match(gene_order, rownames(mat)), , drop = FALSE]

# 4️⃣ Scale expression per gene (row z-score)
mat_scaled <- t(scale(t(mat)))

# 5️⃣ Define CTS colors
cts_cols <- c("CTS1"="royalblue", "CTS2"="#B2DF8A", "CTS3"="orange")
col_fun <- colorRamp2(c(-4,0,4), c("#4B4B9A","white","#C47E32"))  # purple-white-orange

# 6️⃣ Create column annotation for CTS and header labels
cts_annot <- HeatmapAnnotation(
  CTS = pred$CTS,
  col = list(CTS = cts_cols),
  annotation_name_side = "left",
  show_annotation_name = FALSE,
  simple_anno_size = unit(4, "mm")
)

# Get group counts and breaks
cts_counts <- table(pred$CTS)
group_gaps <- cumsum(cts_counts)
cts_labels <- paste0(names(cts_counts), "\n", "n = ", as.numeric(cts_counts))

# 7️⃣ Plot heatmap
pdf("Figure/cts_yourcohort_healthy_heatmap_NatureStyle.pdf", width = 8.5, height = 4)

Heatmap(
  mat_scaled,
  name = "Scaled expression",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 9),
  column_split = pred$CTS,
  top_annotation = cts_annot,
  column_title = paste(cts_labels, collapse = "     "),
  column_title_gp = gpar(fontsize = 10, fontface = "italic"),
  border = FALSE,
  heatmap_legend_param = list(
    at = c(-4, -2, 0, 2, 4),
    title = "Scaled expression",
    legend_width = unit(3, "cm")
  )
)
dev.off()

cat("✅ Saved Nature-style CTS heatmap to Figure/cts_yourcohort_healthy_heatmap_NatureStyle.pdf\n")


# 2) Load healthy logCPM matrix (the expr_logcpm you already created)
expr_logcpm_hv <- expr_logcpm   # from your code

# 3) Find common genes
common_genes <- intersect(rownames(exp_core_g), rownames(expr_logcpm_hv))

core_mat <- exp_core_g[common_genes, , drop = FALSE]
hv_mat   <- expr_logcpm_hv[common_genes, , drop = FALSE]

# 4) Training labels from sepsis
y_core <- factor(core_samples[colnames(core_mat), "CTS"])

# 5) Train RF (NO COMBAT)
set.seed(123)
rf_model <- randomForest(x = t(core_mat), y = y_core)

# 6) Predict healthy CTS BEFORE ComBat
pred_hv_noCombat <- data.frame(
  Sample = colnames(hv_mat),
  CTS = predict(rf_model, newdata = t(hv_mat))
)

# 7) Save
write.csv(pred_hv_noCombat,
          "Figure/cts_yourcohort_healthy_predictions_noComBat.csv",
          row.names = FALSE)

cat("✅ Saved BEFORE-ComBat Healthy CTS to Figure/cts_yourcohort_healthy_predictions_noComBat.csv\n")



library(tidyverse)
before <- read.csv("Figure/cts_yourcohort_healthy_predictions_noComBat.csv")
after  <- read.csv("Figure/cts_yourcohort_healthy_predictions.csv")

before$CTS_before <- factor(before$CTS, levels = c("1","2","3"))
after$CTS_after   <- factor(after$CTS,  levels = c("1","2","3"))

# merge
df <- merge(before[,c("Sample","CTS_before")],
            after[,c("Sample","CTS_after")],
            by = "Sample")

library(dplyr)
library(tidyr)

# df has: Sample, CTS_before, CTS_after
long_df <- df %>%
  tidyr::pivot_longer(
    cols = c(CTS_before, CTS_after),
    names_to = "State",
    values_to = "CTS"
  )

prop_df <- long_df %>%
  dplyr::count(State, CTS, name = "N") %>%   # safer than summarise + n()
  group_by(State) %>%
  mutate(Percent = 100 * N / sum(N))


## ============================================================
## 2. Figure 1 — Stacked bar plot (Before vs After)
## ============================================================

# reorder factor for plotting
prop_df$CTS <- factor(prop_df$CTS, levels = c("1","2","3"))

# Set the desired order: AFTER first, BEFORE second
prop_df$State <- factor(
  prop_df$State,
  levels = c("CTS_after", "CTS_before"),
  labels = c("CTS_ComBat", "CTS_nonCombat")
)

# Colors
cts_colors <- c("1" = "royalblue", "2" = "#B2DF8A", "3" = "orange")

p1 <- ggplot(prop_df, aes(x = State, y = Percent, fill = CTS)) +
  geom_col(width = 0.6, color = "white", size = 0.6) +
  # ❌ removed geom_text (no white percentage labels)
  scale_fill_manual(values = cts_colors) +
  labs(
    x = "",
    y = "Percentage (%)",
    title = "Non-infectious subtype proportions"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0),
    legend.position = "right"
  )

# save figure
dir.create("Figure", showWarnings = FALSE)
ggsave("Figure/CTS_Healthy_Proportion_BeforeAfter_ComBat.svg",
       p1, width = 6, height = 5)

cat("✅ Saved: Figure/CTS_Healthy_Proportion_BeforeAfter_ComBatsvg\n")








## ================================================================
## 0. Load packages
## ================================================================
library(tidyverse)
library(limma)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(msigdbr)
library(fgsea)
library(ggtext)

## ================================================================
## 1. Prepare input: expression + CTS labels for healthy cohort
## ================================================================
# expr_logcpm_hv   -> your logCPM matrix for HV
# res_healthy      -> output from run_subtype_classifier()

CTS_predictions_hv <- res_healthy$predictions
expr_logcpm_df_hv  <- expr_logcpm_hv      # rename to match CAP script
expr <- as.matrix(expr_logcpm_df_hv)

# Keep only HV samples
expr <- expr[, CTS_predictions_hv$Sample, drop = FALSE]

# CTS factor
cts_df <- CTS_predictions_hv %>%
  mutate(CTS = factor(CTS, levels = c("1","2","3"),
                      labels = c("CTS1","CTS2","CTS3")))

# Ensure gene IDs are without version numbers
rownames(expr) <- sub("\\..*$","", rownames(expr))
storage.mode(expr) <- "double"


## ================================================================
## 2. Ensembl → Entrez mapping
## ================================================================
entrez_map <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys       = unique(rownames(expr)),
  keytype    = "ENSEMBL",
  column     = "ENTREZID",
  multiVals  = "first"
)


## ================================================================
## 3. limma: CTS1/2/3 vs others (same as CAP)
## ================================================================
design <- model.matrix(~0 + CTS, data = cts_df)
colnames(design) <- c("CTS1","CTS2","CTS3")

contr <- makeContrasts(
  CTS1_vs_others = CTS1 - (CTS2 + CTS3)/2,
  CTS2_vs_others = CTS2 - (CTS1 + CTS3)/2,
  CTS3_vs_others = CTS3 - (CTS1 + CTS2)/2,
  levels = design
)

fit  <- lmFit(expr, design)
fit2 <- contrasts.fit(fit, contr)
fit3 <- eBayes(fit2)
tmat <- fit3$t    # t-statistic matrix


## ================================================================
## 4. Hallmark gene sets
## ================================================================
hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, entrez_gene)

hallmark_sets <- split(hallmark$entrez_gene, hallmark$gs_name)


## ================================================================
## 5. GSEA geneList function (same as CAP)
## ================================================================
make_geneList <- function(t_vec, ens_ids, entrez_map) {
  tibble(ENSEMBL = ens_ids, t = as.numeric(t_vec)) %>%
    mutate(ENTREZ = entrez_map[ENSEMBL]) %>%
    filter(!is.na(ENTREZ)) %>%
    group_by(ENTREZ) %>%
    slice_max(order_by = abs(t), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    { setNames(.$t, .$ENTREZ) } %>%
    sort(decreasing = TRUE)
}

## ================================================================
## 6. Run GSEA for CTS1/2/3 vs others
## ================================================================
set.seed(1)
gsea_all_hv <- lapply(colnames(tmat), function(cn) {
  gl <- make_geneList(tmat[, cn], rownames(tmat), entrez_map)
  
  fg <- fgsea(
    pathways = hallmark_sets,
    stats    = gl,
    minSize  = 10,
    maxSize  = 500
  ) %>%
    arrange(padj) %>%
    mutate(contrast = cn)
  
  fg
}) %>% bind_rows()

gsea_all_hv <- gsea_all_hv %>% dplyr::rename(gs_name = pathway)


## ================================================================
## 7. Keep best pathway per Hallmark (padj smallest; |NES| to break ties)
## ================================================================
df_best <- gsea_all_hv %>%
  dplyr::mutate(
    pathway      = str_remove(gs_name, "^HALLMARK_") %>% str_replace_all("_"," "),
    contrast     = recode(contrast,
                          CTS1_vs_others = "CTS1",
                          CTS2_vs_others = "CTS2",
                          CTS3_vs_others = "CTS3"),
    FDRneglog10  = -log10(padj)
  ) %>%
  dplyr::filter(!is.na(padj), padj < 0.05) %>%
  group_by(pathway) %>%
  dplyr::arrange(padj, desc(abs(NES)), .by_group = TRUE) %>%
  dplyr::slice(1) %>%
  ungroup()


## ================================================================
## 8. Order & label (same method as CAP)
## ================================================================
df_best2 <- df_best %>%
  dplyr::mutate(
    contrast      = factor(contrast, levels = c("CTS1","CTS2","CTS3")),
    pathway_label = stringr::str_to_sentence(pathway)
  ) %>%
  dplyr::arrange(contrast, desc(FDRneglog10)) %>%
  dplyr::mutate(pathway_ord = factor(pathway_label, levels = rev(pathway_label)))


## ================================================================
## 9. Color palette (same as CAP)
## ================================================================
cts_cols <- c(
  "CTS1"="#355E9A",
  "CTS2"="#9ACD66",
  "CTS3"="#F5A623"
)


## ================================================================
## 10. FINAL FIGURE — HEALTHY GSEA
## ================================================================
p_hv <- ggplot(df_best2, aes(x = FDRneglog10, y = pathway_ord, fill = contrast)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.25) +
  scale_fill_manual(values = cts_cols, name = NULL) +
  labs(x = expression("FDR q ("*-log[10]*")"), y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    axis.text.y  = element_text(size = 10),
    axis.text.x  = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.title = element_text(hjust = 0)
  )

print(p_hv)

ggsave("Figure/CTS_Healthy_GSEA_final.svg",
       p_hv, width = 6.5, height = 7, dpi = 300)

ggsave("Figure/CTS_Healthy_GSEA_final.pdf",
       p_hv, width = 6.5, height = 7)

