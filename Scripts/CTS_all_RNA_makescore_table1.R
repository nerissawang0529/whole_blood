rm(list = ls())

## ===============================
## 1) Install dependencies cleanly
## ===============================
# Fresh start recommended: Session -> Restart R before running.

# Use binaries when possible (Mac-friendly)

## install.packages("pacman")
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
               tableone, kmed, survival, survminer, skimr, lattice, reshape2, grDevices, ggrepel)

options(pkgType = "binary")

# Managers
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("remotes", quietly = TRUE))     install.packages("remotes")

# Make sure your Bioconductor is aligned to your R (R 4.4.x -> Bioc 3.20)
BiocManager::install(version = "3.20", ask = FALSE)

# Install required deps (from Bioconductor/CRAN)
BiocManager::install(c("mixOmics", "sva"), ask = FALSE, update = FALSE)

# Fallback for mixOmics if Bioc binary is unavailable on your machine
if (!requireNamespace("mixOmics", quietly = TRUE)) {
  remotes::install_github("mixOmicsTeam/mixOmics", upgrade = "never")
}

## ==========================
## 2) Install the CTS package
## ==========================
remotes::install_github(
  "bpscicluna/ConsensusTranscriptomicSubtype",
  upgrade = "never",
  dependencies = TRUE,
  build_vignettes = FALSE
)

library(ConsensusTranscriptomicSubtype)


# --- 1) Load the package + your uploaded example objects ---
if (!requireNamespace("ConsensusTranscriptomicSubtype", quietly = TRUE)) {
  stop("Install the package first with remotes::install_github('bpscicluna/ConsensusTranscriptomicSubtype').")
}
library(ConsensusTranscriptomicSubtype)

#load the example data to try the package####
load("Original_data/exp_core_g.rda")   # provides exp_core_g
load("Original_data/core_samples.rda") # provides core_samples

cat("exp_core_g dims: ", paste(dim(exp_core_g), collapse = " x "), "\n")
cat("core_samples rows: ", nrow(core_samples), "\n")

# --- 2) Make a demo 'external' expression matrix from core (genes match perfectly) ---
set.seed(1)
n_new <- min(20, ncol(exp_core_g))                 # choose up to 20 columns
pick  <- sample(colnames(exp_core_g), n_new)       # sample existing columns

new_expr_demo <- exp_core_g[, pick, drop = FALSE]  # genes x samples
# Add small noise so it's not identical
new_expr_demo <- new_expr_demo + matrix(rnorm(length(new_expr_demo), sd = 0.05),
                                        nrow = nrow(new_expr_demo))
# IMPORTANT: make new sample names unique so merge() won’t create .x/.y
colnames(new_expr_demo) <- paste0("EXT_", colnames(new_expr_demo))

## 3) Output folder + file paths
outdir <- "Original_Data"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

heatmap_file     <- file.path(outdir, "cts_heatmap_demo.pdf")
silhouette_file  <- file.path(outdir, "cts_silhouette_demo.pdf")
pred_file        <- file.path(outdir, "CTS_predictions_demo.csv")
new_expr_file    <- file.path(outdir, "new_expr_demo.csv")
core_combat_file <- file.path(outdir, "exp_core_combat_demo.csv")
new_combat_file  <- file.path(outdir, "new_expr_combat_demo.csv")
session_file     <- file.path(outdir, "sessionInfo_CTS_demo.txt")

## 4) Run the ORIGINAL GitHub classifier
res_demo <- run_subtype_classifier(
  new_expr_data   = new_expr_demo,
  exp_core_g      = exp_core_g,
  core_samples    = core_samples,
  heatmap_file    = heatmap_file,
  silhouette_file = silhouette_file
)

## 5) Save ALL requested outputs (data + figures) into Original_Data

# 5a) Predictions (sample → CTS)
write.csv(res_demo$predictions, pred_file, row.names = FALSE)

# 5b) The demo external matrix we used (with 'gene' column first for readability)
new_expr_out <- cbind(gene = rownames(new_expr_demo), as.data.frame(new_expr_demo))
write.csv(new_expr_out, new_expr_file, row.names = FALSE)

# 5c) ComBat-corrected matrices returned by the function
#     Save with 'gene' column so they’re easy to re-load
exp_core_combat <- res_demo$expression_corrected$core
new_expr_combat <- res_demo$expression_corrected$new

exp_core_combat_out <- cbind(gene = rownames(exp_core_combat), as.data.frame(exp_core_combat))
new_expr_combat_out <- cbind(gene = rownames(new_expr_combat), as.data.frame(new_expr_combat))

write.csv(exp_core_combat_out, core_combat_file, row.names = FALSE)
write.csv(new_expr_combat_out, new_combat_file, row.names = FALSE)

# 5d) Save session info for reproducibility
writeLines(capture.output(sessionInfo()), session_file)

## 6) Console summary
cat("\nSaved to:", normalizePath(outdir), "\n")
cat(" -", basename(heatmap_file), "\n")
cat(" -", basename(silhouette_file), "\n")
cat(" -", basename(pred_file), "\n")
cat(" -", basename(new_expr_file), "\n")
cat(" -", basename(core_combat_file), "\n")
cat(" -", basename(new_combat_file), "\n")
cat(" -", basename(session_file), "\n")

## 7) Quick checks
print(head(res_demo$predictions))
if (!is.null(res_demo$silhouette)) {
  sil <- res_demo$silhouette[, "sil_width"]
  cat("Mean silhouette width (demo):", round(mean(sil), 3), "\n")
}


##adapt to my data####
## --- 1) Start from your 'cnt' data frame (like in your screenshot) ---
# the following information come from the code 'Clean_allRNA_allstudydata'
cnt <- read.csv("Original_data/Combined_counts_Amsterdam_my_cohort_allpatients.csv", check.names = FALSE)


# 1) Take the first column as Ensembl IDs, strip version
ens_raw    <- as.character(cnt[[1]])
ens_nover  <- sub("\\..*$", "", ens_raw)   # remove version suffix

# 2) Remove that first column from cnt
cnt[[1]] <- NULL

# 3) Assign stripped IDs as rownames
rownames(cnt) <- ens_nover

# 4) Convert to numeric matrix
cnt_mat <- as.matrix(cnt)
storage.mode(cnt_mat) <- "numeric"

# 5) Handle duplicates (if any genes collapse after version stripping)
if (any(duplicated(rownames(cnt_mat)))) {
  cnt_mat <- do.call(rbind, lapply(split(cnt_mat, rownames(cnt_mat)), function(df) {
    df <- as.matrix(df)
    df[which.max(rowSums(df)), , drop = FALSE]  # keep row with highest total counts
  }))
}

## Quick sanity check
cat("Genes:", nrow(cnt_mat), "  Samples:", ncol(cnt_mat), "\n")
head(rownames(cnt_mat))
head(cnt_mat[,1:5])

#convert to logCPM for CTS
library(edgeR)
dge <- edgeR::DGEList(counts = cnt_mat)
expr_logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)  # genes x samples
write.csv(as.data.frame(expr_logcpm),
          "Original_data/expr_logcpm.csv",
          row.names = TRUE)   # <-- saves rownames as first column

library(ConsensusTranscriptomicSubtype)

# Load the CTS training objects
load("Original_data/exp_core_g.rda")
load("Original_data/core_samples.rda")

# Ensure overlap with the 18 required CTS genes
missing <- setdiff(rownames(exp_core_g), rownames(expr_logcpm))
if (length(missing)) {
  stop(paste("Missing CTS genes:", paste(missing, collapse=", ")))
}

# Order rows exactly like exp_core_g
new_expr <- expr_logcpm[rownames(exp_core_g), , drop = FALSE]

# Run classifier
dir.create("Results", showWarnings = FALSE)
res <- run_subtype_classifier(
  new_expr_data   = new_expr,
  exp_core_g      = exp_core_g,
  core_samples    = core_samples,
  heatmap_file    = "Results/cts_heatmap.pdf",
  silhouette_file = "Results/cts_silhouette.pdf"
)
# --- EXTRA: clearer cohort heatmap with a thick CTS annotation bar ---
suppressPackageStartupMessages(library(pheatmap)); library(RColorBrewer); library(dplyr)

pred_df <- res$predictions %>% dplyr::arrange(CTS)
mat <- res$expression_corrected$new[, pred_df$Sample, drop = FALSE]

ann <- data.frame(CTS = factor(pred_df$CTS,
                               levels = c("1","2","3"),
                               labels = c("CTS-1","CTS-2","CTS-3")))
rownames(ann) <- pred_df$Sample

ann_colors <- list(CTS = c("CTS-1" = "royalblue",
                           "CTS-2" = "#B2DF8A",
                           "CTS-3" = "orange"))

pdf("Results/cts_heatmap_ANNOTATED.pdf", width = 9, height = 6)
pheatmap::pheatmap(
  mat,
  color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "PuOr")))(50),
  cluster_rows = TRUE, cluster_cols = FALSE,
  show_rownames = TRUE, show_colnames = FALSE,
  scale = "row", treeheight_row = 0, border_color = NA,
  annotation_col = ann,
  annotation_colors = ann_colors,
  annotation_legend = TRUE,
  annotation_names_col = TRUE
)
dev.off()

# Save predictions
write.csv(res$predictions, "Original_data/CTS_predictions.csv", row.names = FALSE)
head(res$predictions)

# --- Assumes you already have:
#   - new_expr  : your log-CPM matrix (genes x samples) ordered like rownames(exp_core_g)
#   - exp_core_g, core_samples loaded
#   - library(ConsensusTranscriptomicSubtype) loaded

make_cts_outputs <- function(new_expr, exp_core_g, core_samples, outdir = "Original_Data",
                             prefix = "cts") {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  heatmap_file    <- file.path(outdir, paste0(prefix, "_heatmap.pdf"))
  silhouette_file <- file.path(outdir, paste0(prefix, "_silhouette.pdf"))
  pred_file       <- file.path(outdir, paste0(prefix, "_predictions.csv"))
  sil_csv         <- file.path(outdir, paste0(prefix, "_silhouette_widths.csv"))
  
  # Run the classifier (this call *creates* the PDFs)
  res <- ConsensusTranscriptomicSubtype::run_subtype_classifier(
    new_expr_data   = new_expr,
    exp_core_g      = exp_core_g,
    core_samples    = core_samples,
    heatmap_file    = heatmap_file,
    silhouette_file = silhouette_file
  )
  
  # Save predictions
  utils::write.csv(res$predictions, pred_file, row.names = FALSE)
  
  # Save silhouette widths (if available)
  if (!is.null(res$silhouette)) {
    sil_df <- as.data.frame(res$silhouette)
    # keep only sample id + sil_width if present
    keep <- intersect(colnames(sil_df), c("sample", "sil_width", "cluster"))
    if (length(keep)) sil_df <- sil_df[, keep, drop = FALSE]
    utils::write.csv(sil_df, sil_csv, row.names = FALSE)
  }
  
  # Console summary
  cat("\nCTS outputs saved to:", normalizePath(outdir), "\n")
  cat(" -", basename(heatmap_file), "\n")
  cat(" -", basename(silhouette_file), "\n")
  cat(" -", basename(pred_file), "\n")
  if (!is.null(res$silhouette)) cat(" -", basename(sil_csv), "\n")
  
  # Quick overview
  cat("\nCounts per CTS:\n")
  print(table(res$predictions$CTS), quote = FALSE)
  if (!is.null(res$silhouette)) {
    mw <- mean(res$silhouette[, "sil_width"])
    cat("Mean silhouette width:", round(mw, 3), "\n")
  }
  
  invisible(res)
}

# ===== Run for YOUR cohort and write into Original_Data =====
res_yours <- make_cts_outputs(
  new_expr    = new_expr,       # your logCPM matrix aligned to rownames(exp_core_g)
  exp_core_g  = exp_core_g,
  core_samples= core_samples,
  outdir      = "Original_Data",
  prefix      = "cts_yourcohort"
)


# read predictions (or use res_yours$predictions directly)
pred <- read.csv("Original_Data/cts_yourcohort_predictions.csv", check.names = FALSE)

# counts in numeric order 1-2-3
pred$CTS <- factor(pred$CTS, levels = c("1","2","3"))
df <- as.data.frame(table(pred$CTS))
names(df) <- c("CTS","N")

# bar chart with labels
library(ggplot2)
p <- ggplot(df, aes(x = CTS, y = N)) +
  geom_col() +
  geom_text(aes(label = N), vjust = -0.4) +
  labs(title = "CTS subtype counts (your cohort)",
       x = "CTS subtype", y = "Number of samples") +
  theme_classic()

# show on screen
print(p)

# save to PDF
ggsave("Original_Data/cts_yourcohort_subtype_counts_ggplot.pdf", p, width = 5, height = 4)



#combine the CTS score in the table 1 dataframe
dfwb_sel <- read.csv("Original_Data/qns_ID.csv", check.names = FALSE)
# the dataframe 'qns_ID' is from code 'Clean_allRNA_allstudydata' and dataframe name is 'Elder_cap'
library(dplyr)


merged_df <- pred %>%
  dplyr::left_join(dfwb_sel %>% dplyr::select(Sample, ID), by = "Sample")



#table one for CTS
studydata <-read.csv("Original_Data/studydata_matched_allpatients_tcs.csv", check.names = FALSE)
# Keep only rows in studydata that appear in merged_df$ID and add CTS (base R)
cts <- unique(merged_df[c("ID","CTS")]); names(cts)[1] <- "EB_id"
studydata_with_cts <- merge(transform(studydata, EB_id = as.character(EB_id)),
                            transform(cts,       EB_id = as.character(EB_id)),
                            by = "EB_id", all = FALSE)

destination_folder <- "Original_data/" 
export_file_name <- "studydata_with_cts.csv" 
write.csv(studydata_with_cts, file.path(destination_folder, export_file_name), row.names = FALSE)
cat("File exported to:", file.path(destination_folder, export_file_name), "\n")

## table1
## analyse errors 
allvars <- c("group","flow_Data", "spectral", "inclusion_hospital", "sampling_time",		"age_yrs",	"gender",	"ethnic_group#White",	
             "BMI",	"symptom_days",
             "CCI",	"hypertension",	"cpd",	"COPD",	"diabetes",	
             "ccd","ckd",	"mneoplasm",	"immune_sup",	"cnd",	
             "qSOFA_score",	"MEWS_score",	"CURB_score",	"PSI_new",
             "crp_1_1",	 "WBC_2_1",	"Neutrophil_unit_1",	"Lymphocyte_1_1",	
             "Creatinine_value_1",
             "Platelets_value_1",	
             "Blood_Urea_Nitrogen_value_1",	
             "Oxygen_therapy_1", 
             "length_of_oxygen", "antibiotic_seven_days", 
             "length_of_stay", 
             "ICU_Medium_Care_admission_1",	 "icu_stay", 
             "Non_invasive_ventilation_1",	"Invasive_ventilation_1",	"lenght_of_intubation",
             "bacteremia", "pathogen_cultured",
             "hospdeath",	"mortality_d30", "mortality_d90","CTS","ttcs_halm_372_days")

catvars <- c("group","flow_Data", "spectral", "inclusion_hospital", "gender",	"ethnic_group#White",	
             "hypertension",	"cpd",	"COPD",	"diabetes",	
             "ccd","ckd",	"mneoplasm",	"immune_sup",	"cnd",	
             "Oxygen_therapy_1", 
             "antibiotic_seven_days", 
             "ICU_Medium_Care_admission_1",	 
             "Non_invasive_ventilation_1",	"Invasive_ventilation_1",	
             "bacteremia",
             "hospdeath",	"mortality_d30", "mortality_d90", "pathogen_cultured","CTS")

nonnormal <- c("sampling_time",		"age_yrs", "BMI", "symptom_days",  "CCI", "qSOFA_score",	"MEWS_score",	"CURB_score",	"PSI_new",
               "crp_1_1",	 "WBC_2_1",	"Neutrophil_unit_1",	"Lymphocyte_1_1",	
               "Creatinine_value_1",
               "Platelets_value_1",	
               "Blood_Urea_Nitrogen_value_1",	"length_of_oxygen","length_of_stay" ,"icu_stay", "lenght_of_intubation","ttcs_halm_372_days")

tab1     <- CreateTableOne(
  vars        = allvars, 
  data        = studydata_with_cts, 
  strata = "CTS",
  factorVars  = catvars,
  test        = TRUE)

print(tab1, nonnormal = nonnormal, quote = FALSE, noSpaces = TRUE, smd = T, missing = T)

# Convert the CreateTableOne object to a data frame
tab1_df <- as.data.frame(print(tab1, nonnormal = nonnormal, quote = FALSE, noSpaces = TRUE, smd = TRUE, missing = TRUE))
tab1_df$tablename <- rownames(tab1_df) 
tab1_df <- tab1_df[, c("tablename", setdiff(names(tab1_df), "tablename"))]

# Write the data frame to a CSV file
#export clinical_marker_unique data
destination_folder <- "Original_data/" 
export_file_name <- "table1_cts.csv" 
write.csv(tab1_df, file.path(destination_folder, export_file_name), row.names = FALSE)
cat("File exported to:", file.path(destination_folder, export_file_name), "\n")



##heatmap####
# --- build matrix + annotation as you already do ---
pred_df <- res$predictions %>% dplyr::arrange(CTS)
library(org.Hs.eg.db)
library(AnnotationDbi)

# --- Convert Ensembl IDs to gene symbols ---
ensembl_ids <- rownames(res$expression_corrected$new)
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = sub("\\..*$", "", ensembl_ids),  # 去掉版本号
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")

# 替换 rownames 中的 Ensembl 为基因名
rownames(res$expression_corrected$new) <- ifelse(is.na(gene_symbols),
                                                 ensembl_ids,
                                                 gene_symbols)

# 确保后续 heatmap 用的是 gene symbol
mat <- res$expression_corrected$new[, pred_df$Sample, drop = FALSE]

ann <- data.frame(
  CTS = factor(pred_df$CTS, levels = c("1","2","3"),
               labels = c("CTS-1","CTS-2","CTS-3"))
)
rownames(ann) <- pred_df$Sample

ann_colors <- list(CTS = c("CTS-1" = "royalblue", "CTS-2" = "#B2DF8A", "CTS-3" = "orange"))

# --- palette + symmetric breaks around 0 (z-scores) ---
library(RColorBrewer)
pal     <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(101)  # blue-white-red
brks    <- seq(-4, 4, length.out = length(pal))                # clamp at [-4, 4]

pdf("Original_Data/cts_yourcohort_heatmap_ANNOTATED.pdf", width = 9, height = 6)
# --- Reorder genes to match the paper figure ---
## ---- lock gene order just before plotting ----
desired_order <- c(
  "ACER3","SERPINB1","HK3","TDRD9","NLRC4","PGD",
  "UBE2H","METTL9","STOM","SNX3","GADD45A","BTN3A3",
  "BPGM","CA1","SLC4A1","EPB42","FECH","GLRX5"
)

# 只保留出现的基因，并严格按 desired_order 排（不会产生 NA）
idx <- match(desired_order, rownames(mat))
idx <- idx[!is.na(idx)]
mat <- mat[idx, , drop = FALSE]

# （可选）检查一下当前行顺序
cat("Row order:\n"); print(rownames(mat))
pdf("Results/cts_heatmap_ANNOTATED_ordered.pdf", width = 9, height = 6)
pheatmap::pheatmap(
  mat,
  color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "PuOr")))(50),
  cluster_rows = FALSE,     # ← 关键：不要再聚类行
  cluster_cols = FALSE,
  show_rownames = TRUE, show_colnames = FALSE,
  scale = "row", treeheight_row = 0, border_color = NA,
  annotation_col = ann,
  annotation_colors = ann_colors,
  annotation_legend = TRUE,
  annotation_names_col = TRUE
)
dev.off()


pheatmap::pheatmap(
  mat,
  scale = "row",                                  # z-score per gene (matches pkg heatmap)
  color = pal, breaks = brks,                     # mimic the first figure’s look
  cluster_rows = TRUE, cluster_cols = FALSE,
  show_rownames = TRUE, show_colnames = FALSE,
  border_color = NA, treeheight_row = 0,
  annotation_col = ann, annotation_colors = ann_colors,
  annotation_legend = TRUE, annotation_names_col = TRUE
)
dev.off()
##
suppressPackageStartupMessages({
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
})

# 
pred_df <- res$predictions %>% dplyr::arrange(CTS)
mat     <- res$expression_corrected$new[, pred_df$Sample, drop = FALSE]

ann <- data.frame(
  CTS = factor(pred_df$CTS, levels = c("1","2","3"),
               labels = c("CTS-1","CTS-2","CTS-3"))
)
rownames(ann) <- pred_df$Sample

ann_colors <- list(CTS = c("CTS-1" = "royalblue",
                           "CTS-2" = "#B2DF8A",
                           "CTS-3" = "orange"))

# --- palette + symmetric breaks around 0 (z-scores) ---
pal  <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(101)
brks <- seq(-4, 4, length.out = length(pal))

# 
if (!requireNamespace("svglite", quietly = TRUE)) install.packages("svglite")
svglite::svglite("Original_Data/cts_yourcohort_heatmap_ANNOTATED.svg",
                 width = 9, height = 6)

pheatmap::pheatmap(
  mat,
  scale = "row",
  color = pal, breaks = brks,
  cluster_rows = TRUE, cluster_cols = FALSE,
  show_rownames = TRUE, show_colnames = FALSE,
  border_color = NA, treeheight_row = 0,
  annotation_col = ann, annotation_colors = ann_colors,
  annotation_legend = TRUE, annotation_names_col = TRUE
)
dev.off()


## 直接跑分类器（结果保存在 Results 目录）
dir.create("Original_Data", showWarnings = FALSE)

res <- run_subtype_classifier(
  new_expr_data   = new_expr,
  exp_core_g      = exp_core_g,
  core_samples    = core_samples,
  heatmap_file    = "Original_Data/cts_heatmap.pdf",   # 保持不改（你没要求改）
  silhouette_file = "Original_Data/cts_silhouette.svg" # ← 这里改成 .svg
)

make_cts_outputs <- function(new_expr, exp_core_g, core_samples,
                             outdir = "Original_Data", prefix = "cts") {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  heatmap_file    <- file.path(outdir, paste0(prefix, "_heatmap.pdf"))      # 不改
  silhouette_file <- file.path(outdir, paste0(prefix, "_silhouette.svg"))   # ← 改为 .svg
  pred_file       <- file.path(outdir, paste0(prefix, "_predictions.csv"))
  sil_csv         <- file.path(outdir, paste0(prefix, "_silhouette_widths.csv"))
  
  res <- ConsensusTranscriptomicSubtype::run_subtype_classifier(
    new_expr_data   = new_expr,
    exp_core_g      = exp_core_g,
    core_samples    = core_samples,
    heatmap_file    = heatmap_file,
    silhouette_file = silhouette_file
  )
  
  utils::write.csv(res$predictions, pred_file, row.names = FALSE)
  
  if (!is.null(res$silhouette)) {
    sil_df <- as.data.frame(res$silhouette)
    keep <- intersect(colnames(sil_df), c("sample", "sil_width", "cluster"))
    if (length(keep)) sil_df <- sil_df[, keep, drop = FALSE]
    utils::write.csv(sil_df, sil_csv, row.names = FALSE)
  }
  
  cat("\nCTS outputs saved to:", normalizePath(outdir), "\n")
  cat(" -", basename(heatmap_file), "\n")
  cat(" -", basename(silhouette_file), "\n")
  cat(" -", basename(pred_file), "\n")
  if (!is.null(res$silhouette)) cat(" -", basename(sil_csv), "\n")
  
  cat("\nCounts per CTS:\n")
  print(table(res$predictions$CTS), quote = FALSE)
  if (!is.null(res$silhouette)) {
    mw <- mean(res$silhouette[, "sil_width"])
    cat("Mean silhouette width:", round(mw, 3), "\n")
  }
  invisible(res)
}
