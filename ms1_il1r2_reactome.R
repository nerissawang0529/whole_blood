############################################################
## Reactome GSEA: MS1 vs IL1R2+ SIDE-BY-SIDE BAR PLOT
## - Mother pathways: biologically predefined backbone (reactome_gl.txt)
## - Child pathways: UNION rule (GLOBAL BH FDR < 0.05 in either analysis)
## - Main figure uses GLOBAL BH only (ReactomePA p.adjust); NO subset-BH
## - Color rule: p.adjust < 0.05 -> NES color (red/blue); otherwise grey
## - Mother pathways ALWAYS display a bar (even if not significant)
## - MS1 left, IL1R2 right
##
## Targeted updates:
##  1) Circularity control: mapIds(..., multiVals="list") + unlist -> remove ALL mapped ENSEMBL
##  2) fgsea warnings: break ties in gene-level stats (tiny deterministic jitter) + eps=0
##  3) Global style locked: Arial + sizes + line widths + margins
############################################################

############################################################
## 0) Clean start
############################################################
rm(list = ls())

############################################################
## 1) Packages
############################################################
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

pkgs_cran <- c(
  "dplyr","readr","tibble","stringr","tidyr","purrr",
  "ggplot2","scales","ggtext","patchwork"
)
pkgs_bioc <- c(
  "DESeq2","ReactomePA","org.Hs.eg.db","AnnotationDbi"
)

pacman::p_load(char = c(pkgs_cran, pkgs_bioc))

############################################################
## Global style (LOCKED)
############################################################
FONT_FAMILY <- "Arial"  # Font family used for all text
BASE_SIZE   <- 14       # Base font size (pt)
TITLE_SIZE  <- 14       # Plot title size (pt)
TAG_SIZE    <- 16       # Panel tag size (pt) - NOT bold
AXIS_TITLE  <- 13       # Axis title size (pt)
AXIS_TEXT   <- 11       # Axis tick label size (pt)

LINE_AXIS   <- 0.6      # Axis line width
LINE_TICKS  <- 0.6      # Tick line width
BIN_EDGE_W  <- 0.2      # Histogram bin edge width
MARGINS_MM  <- c(5.5, 5.5, 5.5, 5.5)  # Plot margins (t, r, b, l) in mm

############################################################
## 2) Inputs (files)
############################################################
counts_file    <- "Original_Data/Combined_counts_Amsterdam_my_cohort_allpatients.csv"
pheno_file     <- "Original_Data/studydata_with_ms1_il1r2.csv"
ms1_sig_file   <- "Original_Data/MS1_septiStates_150_Reyes_2020.txt"
il1r2_sig_file <- "Original_data/CIBERSORTx_Job16_ref_matrix_equal_inferred_phenoclasses.CIBERSORTx_Job16_ref_matrix_equal_inferred_refsample.bm.K999.txt"
edges_file     <- "Original_Data/reactome_gl.txt"

stopifnot(
  file.exists(counts_file),
  file.exists(pheno_file),
  file.exists(ms1_sig_file),
  file.exists(il1r2_sig_file),
  file.exists(edges_file)
)

############################################################
## 3) User-locked settings
############################################################
# Comparisons (locked)
ms1_levels   <- c("MS1_low", "MS1_high")                 # MS1_high vs MS1_low
il1r2_levels <- c("IL1R2_zero", "IL1R2_positive")        # IL1R2_positive vs IL1R2_zero

# Mother pathways (biological definition)
mother_keep <- c(
  "Innate Immune System",
  "Adaptive Immune System",
  "Cytokine Signaling in Immune system",
  "Hemostasis"
)

# GLOBAL BH threshold
fdr_cutoff_union <- 0.05

# Optional fallback (only when a mother would be empty in DISPLAY)
# (kept, but now "significance" still comes from global p.adjust; fallback is only for display)
fallback_n_max     <- 2
fallback_nominal_p <- 0.05

# Plot settings
bar_width <- 0.72

############################################################
## 4) Helper: load counts as integer matrix (FIXED)
############################################################
load_counts_int <- function(file) {
  cnt <- read.csv(file, check.names = FALSE)
  
  # First column = Ensembl IDs (possibly with version)
  ens_raw   <- as.character(cnt[[1]])
  ens_clean <- sub("\\..*$", "", ens_raw)
  
  cnt[[1]] <- NULL
  rownames(cnt) <- ens_clean
  
  # Force matrix
  cnt <- as.matrix(cnt)
  
  # Checks
  storage.mode(cnt) <- "numeric"
  if (any(cnt < 0, na.rm = TRUE)) stop("Counts contain negative values.")
  if (any(is.na(cnt))) stop("Counts contain NA values.")
  
  # DESeq2 requires integer counts
  cnt_int <- round(cnt)
  if (!all(cnt_int == cnt)) stop("Counts are not integer-like. DESeq2 requires raw integer counts.")
  storage.mode(cnt_int) <- "integer"
  
  cnt_int
}

############################################################
## 5) Helper: align phenotype to counts
############################################################
align_meta_counts <- function(cnt_int, studydata, group_col, levels_vec) {
  studydata <- studydata %>% dplyr::mutate(Patient_ID = as.character(Patient_ID))
  
  required_cols <- c("Patient_ID", group_col)
  missing_cols  <- setdiff(required_cols, colnames(studydata))
  if (length(missing_cols) > 0) {
    stop(paste0("Phenotype file missing columns: ", paste(missing_cols, collapse = ", ")))
  }
  
  missing_in_counts <- setdiff(studydata$Patient_ID, colnames(cnt_int))
  if (length(missing_in_counts) > 0) {
    stop(paste0(
      "These Patient_ID are in phenotype but not in counts: ",
      paste(head(missing_in_counts, 30), collapse = ", "),
      ifelse(length(missing_in_counts) > 30, " ...", "")
    ))
  }
  
  cnt_int <- cnt_int[, studydata$Patient_ID, drop = FALSE]
  stopifnot(all(colnames(cnt_int) == studydata$Patient_ID))
  
  meta <- studydata %>%
    dplyr::transmute(
      Patient_ID = Patient_ID,
      group = factor(.data[[group_col]], levels = levels_vec)
    )
  
  if (any(is.na(meta$group))) {
    stop(paste0(group_col, " has unexpected values. Expected levels: ", paste(levels_vec, collapse = ", ")))
  }
  
  list(cnt_int = cnt_int, meta = meta)
}

############################################################
## 6) Helper: circularity control (remove top1000 signature genes)
##    Targeted update: mapIds multiVals = "list" + unlist
############################################################
remove_signature_top1000 <- function(cnt_int, sig_file, mode = c("MS1", "IL1R2")) {
  mode <- match.arg(mode)
  
  if (mode == "MS1") {
    sig <- read.delim(sig_file, stringsAsFactors = FALSE)
    req <- c("gennames", "MS1")
    if (!all(req %in% colnames(sig))) stop("MS1 signature file must contain columns: gennames, MS1")
    
    top_syms <- sig %>%
      dplyr::arrange(dplyr::desc(abs(MS1))) %>%
      dplyr::slice_head(n = 1000) %>%
      dplyr::pull(gennames)
  }
  
  if (mode == "IL1R2") {
    sig <- read.delim(sig_file, stringsAsFactors = FALSE, check.names = FALSE)
    
    gene_col  <- "NAME"
    score_col <- "IL1R2+_immature_neutrophils"
    
    if (!all(c(gene_col, score_col) %in% colnames(sig))) {
      stop(
        paste0(
          "IL1R2 signature file does not contain expected columns.\n",
          "Expected gene_col  = '", gene_col, "'\n",
          "Expected score_col = '", score_col, "'\n\n",
          "Found columns:\n",
          paste(colnames(sig), collapse = ", ")
        )
      )
    }
    
    top_syms <- sig %>%
      dplyr::transmute(
        gene  = as.character(.data[[gene_col]]),
        score = suppressWarnings(as.numeric(.data[[score_col]]))
      ) %>%
      dplyr::filter(!is.na(gene), gene != "", !is.na(score)) %>%
      dplyr::arrange(dplyr::desc(abs(score))) %>%
      dplyr::slice_head(n = 1000) %>%
      dplyr::pull(gene)
  }
  
  # --- targeted change here ---
  top_ens_list <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys      = unique(top_syms),
    keytype   = "SYMBOL",
    column    = "ENSEMBL",
    multiVals = "list"
  )
  top_ens <- unique(na.omit(unlist(top_ens_list)))
  # ---------------------------
  
  cnt_int[!rownames(cnt_int) %in% top_ens, , drop = FALSE]
}

############################################################
## 7) Helper: run DESeq2 + Reactome GSEA and return reactome_df
##    Targeted update: break ties in stats + eps=0 to silence warnings
############################################################
run_deseq2_reactome <- function(cnt_int, meta, contrast_vec) {
  
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = cnt_int,
    colData   = meta %>% dplyr::select(group),
    design    = ~ group
  )
  
  # Low count filtering
  keep <- rowSums(DESeq2::counts(dds) >= 10) >= 3
  dds  <- dds[keep, ]
  
  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::results(dds, contrast = contrast_vec)
  
  res_df <- as.data.frame(res)
  res_df$ENSEMBL <- rownames(res_df)
  res_df <- res_df %>% dplyr::filter(!is.na(stat))
  
  res_df$entrez <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys      = res_df$ENSEMBL,
    keytype   = "ENSEMBL",
    column    = "ENTREZID",
    multiVals = "first"
  )
  
  gene_df <- res_df %>%
    dplyr::filter(!is.na(entrez)) %>%
    dplyr::group_by(entrez) %>%
    dplyr::slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
  
  gene_list <- gene_df$stat
  names(gene_list) <- gene_df$entrez
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  # --- targeted change: break ties deterministically (avoids "ties in preranked stats") ---
  if (any(duplicated(gene_list))) {
    set.seed(11)  # deterministic
    gene_list <- gene_list + runif(length(gene_list), min = -1e-7, max = 1e-7)
    gene_list <- sort(gene_list, decreasing = TRUE)
  }
  # --------------------------------------------------------------------------------------
  
  set.seed(11)
  gsea <- ReactomePA::gsePathway(
    geneList      = gene_list,
    organism      = "human",
    pvalueCutoff  = 1,
    pAdjustMethod = "BH",   # GLOBAL BH (this is what we will use)
    minGSSize     = 10,
    maxGSSize     = 10000000,
    verbose       = FALSE,
    nPermSimple   = 10000,
    eps           = 0       # targeted change: better estimation for extremely small p-values
  )
  
  as.data.frame(gsea)
}

############################################################
## 8) Prepare Reactome backbone (Mother / Children)
############################################################
reactome_edges <- read.delim(edges_file) %>%
  dplyr::mutate(
    Mother   = stringr::str_trim(Mother),
    Children = stringr::str_trim(Children)
  )

if (!all(c("Mother","Children") %in% colnames(reactome_edges))) {
  stop("reactome_gl.txt must contain columns: Mother and Children")
}

reactome_edges <- reactome_edges %>%
  dplyr::filter(Mother %in% mother_keep)

# LOCK mother order exactly as mother_keep
mother_levels <- mother_keep[mother_keep %in% reactome_edges$Mother]

############################################################
## 9) Load phenotype once
############################################################
studydata <- readr::read_csv(pheno_file, show_col_types = FALSE) %>%
  dplyr::mutate(Patient_ID = as.character(Patient_ID))

############################################################
## 10) Run MS1 Reactome (with circularity control)
############################################################
cnt0 <- load_counts_int(counts_file)

aligned_ms1 <- align_meta_counts(cnt0, studydata, group_col = "MS1_group", levels_vec = ms1_levels)
cnt_ms1     <- aligned_ms1$cnt_int
meta_ms1    <- aligned_ms1$meta

cnt_ms1_cc <- remove_signature_top1000(cnt_ms1, ms1_sig_file, mode = "MS1")

reactome_df_ms1 <- run_deseq2_reactome(
  cnt_int      = cnt_ms1_cc,
  meta         = meta_ms1,
  contrast_vec = c("group", "MS1_high", "MS1_low")
)

############################################################
## 11) Run IL1R2 Reactome (with circularity control)
############################################################
aligned_il1 <- align_meta_counts(cnt0, studydata, group_col = "IL1R2_group", levels_vec = il1r2_levels)
cnt_il1     <- aligned_il1$cnt_int
meta_il1    <- aligned_il1$meta

cnt_il1_cc <- remove_signature_top1000(cnt_il1, il1r2_sig_file, mode = "IL1R2")

reactome_df_il1r2 <- run_deseq2_reactome(
  cnt_int      = cnt_il1_cc,
  meta         = meta_il1,
  contrast_vec = c("group", "IL1R2_positive", "IL1R2_zero")
)

############################################################
## 12) Prepare GSEA tables (global p.adjust + pvalue + NES)
############################################################
prep_gsea_tbl <- function(reactome_df, label) {
  reactome_df %>%
    dplyr::mutate(
      Description = stringr::str_trim(Description),
      NES         = as.numeric(NES),
      p.adjust    = as.numeric(p.adjust),   # GLOBAL BH from ReactomePA
      pvalue      = as.numeric(pvalue)
    ) %>%
    dplyr::select(Description, NES, pvalue, p.adjust) %>%
    dplyr::mutate(analysis = label)
}

gsea_ms1 <- prep_gsea_tbl(reactome_df_ms1,   "MS1")
gsea_il1 <- prep_gsea_tbl(reactome_df_il1r2, "IL1R2")

############################################################
## 13) Child pathway selection: UNION rule using GLOBAL BH (p.adjust)
############################################################
all_children <- reactome_edges %>% dplyr::pull(Children) %>% unique()

sig_children_ms1 <- gsea_ms1 %>%
  dplyr::filter(Description %in% all_children, !is.na(p.adjust), p.adjust < fdr_cutoff_union) %>%
  dplyr::pull(Description)

sig_children_il1 <- gsea_il1 %>%
  dplyr::filter(Description %in% all_children, !is.na(p.adjust), p.adjust < fdr_cutoff_union) %>%
  dplyr::pull(Description)

master_children <- sort(unique(c(sig_children_ms1, sig_children_il1)))

############################################################
## 14) Optional fallback: if a Mother would be empty, add up to 2 children
##     with nominal p < 0.05 in either analysis (pvalue),
##     (Coloring still follows GLOBAL BH p.adjust, i.e., may remain grey.)
############################################################
fallback_children <- c()

for (m in mother_levels) {
  children_m <- reactome_edges %>% dplyr::filter(Mother == m) %>% dplyr::pull(Children) %>% unique()
  shown_m    <- intersect(children_m, master_children)
  
  if (length(shown_m) == 0) {
    
    cand_ms1 <- gsea_ms1 %>%
      dplyr::filter(Description %in% children_m, !is.na(pvalue), pvalue < fallback_nominal_p) %>%
      dplyr::arrange(pvalue) %>%
      dplyr::pull(Description)
    
    cand_il1 <- gsea_il1 %>%
      dplyr::filter(Description %in% children_m, !is.na(pvalue), pvalue < fallback_nominal_p) %>%
      dplyr::arrange(pvalue) %>%
      dplyr::pull(Description)
    
    cand_union <- unique(c(cand_ms1, cand_il1))
    
    if (length(cand_union) > 0) {
      fallback_children <- c(fallback_children, head(cand_union, fallback_n_max))
    }
  }
}

fallback_children <- sort(unique(fallback_children))
master_children2  <- sort(unique(c(master_children, fallback_children)))

############################################################
## 15) Build backbone rows (Mother headers + selected Children)
############################################################
df_mothers <- tibble::tibble(
  parent_cat  = mother_levels,
  Description = mother_levels,
  is_mother   = TRUE
)

df_children <- reactome_edges %>%
  dplyr::filter(Children %in% master_children2) %>%
  dplyr::transmute(
    parent_cat  = Mother,
    Description = Children,
    is_mother   = FALSE
  )

df_backbone <- purrr::map_dfr(mother_levels, function(m) {
  dplyr::bind_rows(
    df_mothers  %>% dplyr::filter(parent_cat == m),
    df_children %>% dplyr::filter(parent_cat == m)
  )
})

############################################################
## 16) Join GSEA results (GLOBAL BH) to backbone
##     - Mother ALWAYS shows a bar: if NES is missing, force NES_plot = 0
##     - Significance is based on GLOBAL p.adjust (for BOTH mother and child)
############################################################
join_global <- function(df_backbone, gsea_tbl_one, analysis_label) {
  
  out <- df_backbone %>%
    dplyr::left_join(
      gsea_tbl_one %>% dplyr::select(Description, NES, pvalue, p.adjust),
      by = "Description"
    ) %>%
    dplyr::mutate(
      analysis = analysis_label,
      
      # Mother label bold
      y_label = ifelse(is_mother, paste0("**", parent_cat, "**"), Description),
      
      # Force mother to always have a bar (even if GSEA did not return it)
      NES_plot = dplyr::case_when(
        is_mother & is.na(NES) ~ 0,   # ensure a bar exists
        TRUE ~ NES
      ),
      
      # Global significance for BOTH mother and child
      sig_global = dplyr::case_when(
        !is.na(p.adjust) & p.adjust < 0.05 ~ "FDR < 0.05",
        TRUE ~ "Not significant"
      )
    )
  
  out
}

df_ms1 <- join_global(df_backbone, gsea_ms1, "MS1")
df_il1 <- join_global(df_backbone, gsea_il1, "IL1R2")

df_both <- dplyr::bind_rows(df_ms1, df_il1)

# facet order: MS1 left, IL1R2 right
df_both$analysis <- factor(df_both$analysis, levels = c("MS1", "IL1R2"))

# y order = backbone order, reversed for plotting
y_levels <- rev(unique(df_both$y_label))
df_both$y_label <- factor(df_both$y_label, levels = y_levels)

# Shared NES color range (use NES_plot so mothers with 0 don't break)
nes_range <- range(df_both$NES_plot, na.rm = TRUE)
nes_lim   <- max(abs(nes_range))
if (!is.finite(nes_lim) || nes_lim == 0) nes_lim <- 1

# Fill rule: significant -> NES; otherwise grey (applies to BOTH mother and child)
df_both <- df_both %>%
  dplyr::mutate(
    bar_fill = ifelse(!is.na(NES_plot) & sig_global == "FDR < 0.05", NES_plot, NA_real_)
  )

############################################################
## 17) Plot: side-by-side BAR (GLOBAL BH) with desired theme + LOCKED style
############################################################
p_main_globalBH <- ggplot2::ggplot(df_both, ggplot2::aes(x = NES_plot, y = y_label)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 1.05, color = "black") +
  ggplot2::geom_col(
    ggplot2::aes(fill = bar_fill),
    width = bar_width,
    color = "black",
    linewidth = 0.25,
    na.rm = FALSE
  ) +
  ggplot2::facet_wrap(~analysis, nrow = 1, scales = "free_x") +
  ggplot2::scale_fill_gradient2(
    low  = "#4575b4",
    mid  = "white",
    high = "#d73027",
    limits = c(-nes_lim, nes_lim),
    name = "NES",
    na.value = "grey75"
  ) +
  ggplot2::theme_classic(base_size = BASE_SIZE, base_family = FONT_FAMILY) +
  ggplot2::theme(
    text               = ggplot2::element_text(family = FONT_FAMILY),
    plot.title         = ggplot2::element_text(size = TITLE_SIZE, family = FONT_FAMILY),
    axis.text.y        = ggtext::element_markdown(size = AXIS_TEXT, family = FONT_FAMILY),
    axis.text.x        = ggplot2::element_text(size = AXIS_TEXT, family = FONT_FAMILY),
    axis.title.x       = ggplot2::element_text(size = AXIS_TITLE, family = FONT_FAMILY),
    axis.title.y       = ggplot2::element_text(size = AXIS_TITLE, family = FONT_FAMILY),
    
    axis.ticks.y       = ggplot2::element_blank(),
    axis.ticks.x       = ggplot2::element_line(linewidth = LINE_TICKS, color = "black"),
    axis.line.y        = ggplot2::element_blank(),
    axis.line.x        = ggplot2::element_line(linewidth = LINE_AXIS, color = "black"),
    
    panel.grid         = ggplot2::element_blank(),
    strip.background   = ggplot2::element_blank(),
    strip.text         = ggplot2::element_text(face = "bold", size = BASE_SIZE, family = FONT_FAMILY),
    legend.position    = "right",
    legend.title       = ggplot2::element_text(size = AXIS_TITLE, family = FONT_FAMILY),
    legend.text        = ggplot2::element_text(size = AXIS_TEXT, family = FONT_FAMILY),
    
    axis.title         = ggplot2::element_blank(),
    plot.margin        = ggplot2::margin(
      t = MARGINS_MM[1], r = MARGINS_MM[2], b = MARGINS_MM[3], l = MARGINS_MM[4],
      unit = "mm"
    )
  )
p_main_globalBH <- p_main_globalBH +
  labs(x = "Normalized Enrichment Score (NES)", y = NULL)

print(p_main_globalBH)

############################################################
## 18) OPTIONAL: Save
############################################################

install.packages("remotes")

remotes::install_github("wilkelab/gridtext")
remotes::install_github("wilkelab/ggtext")

library(ggplot2)
library(ggtext)
library(svglite)

methods("element_grob")  

dir.create("Figure", showWarnings = FALSE, recursive = TRUE)

ggsave(
  filename = "Figure/Reactome_MS1_vs_IL1R2_bar_globalBH_main_movegenesignature.svg",
  plot     = p_main_globalBH,
  device   = svglite::svglite,
  width    = 12,
  height   = 8,
  units    = "in"
)

############################################################
## END
############################################################




############################################################
## Reactome GSEA: IL1R2 within MS1_high (after UNION signature removal)
## - Population: MS1_high only
## - Contrast: IL1R2_positive vs IL1R2_zero (reference = IL1R2_zero)
## - Circularity control: remove UNION of top1000 MS1 signature genes
##   and top1000 IL1R2 signature genes BEFORE running DESeq2/GSEA
## - Backbone: biologically predefined Mother/Child map (reactome_gl.txt)
## - Child pathways shown: union rule against GLOBAL BH (p.adjust < 0.05)
## - Mother pathways ALWAYS display a bar (even if missing in GSEA -> NES_plot = 0)
## - Color rule: p.adjust < 0.05 -> NES color (red/blue); otherwise grey
############################################################

#exclude ms1 and il1r2 signatures within MS1_high####
############################################################
############################################################
## Reactome GSEA: IL1R2 within MS1_high (BAR PLOT)
## Circularity control: remove UNION (MS1 top1000 ∪ IL1R2 top1000)
## - mapping rigor: mapIds(..., multiVals="list") + unlist -> remove ALL mapped ENSEMBL
## - fgsea warnings: break ties in gene-level stats (tiny deterministic jitter) + eps=0
## - Global style locked: Arial + sizes + line widths + margins
############################################################

############################################################
## 0) Clean start
############################################################
rm(list = ls())

############################################################
## 1) Packages
############################################################
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

pkgs_cran <- c(
  "dplyr","readr","tibble","stringr","tidyr","purrr",
  "ggplot2","scales","ggtext","svglite"
)
pkgs_bioc <- c(
  "DESeq2","ReactomePA","org.Hs.eg.db","AnnotationDbi"
)

pacman::p_load(char = c(pkgs_cran, pkgs_bioc))

############################################################
## Global style (LOCKED)
############################################################
FONT_FAMILY <- "Arial"
BASE_SIZE   <- 14
TITLE_SIZE  <- 14
TAG_SIZE    <- 16
AXIS_TITLE  <- 13
AXIS_TEXT   <- 11

LINE_AXIS   <- 0.6
LINE_TICKS  <- 0.6
BIN_EDGE_W  <- 0.2
MARGINS_MM  <- c(5.5, 5.5, 5.5, 5.5)

############################################################
## 2) Inputs (files)
############################################################
counts_file    <- "Original_Data/Combined_counts_Amsterdam_my_cohort_allpatients.csv"
pheno_file     <- "Original_Data/studydata_with_ms1_il1r2.csv"
ms1_sig_file   <- "Original_Data/MS1_septiStates_150_Reyes_2020.txt"
il1r2_sig_file <- "Original_data/CIBERSORTx_Job16_ref_matrix_equal_inferred_phenoclasses.CIBERSORTx_Job16_ref_matrix_equal_inferred_refsample.bm.K999.txt"
edges_file     <- "Original_Data/reactome_gl.txt"

stopifnot(
  file.exists(counts_file),
  file.exists(pheno_file),
  file.exists(ms1_sig_file),
  file.exists(il1r2_sig_file),
  file.exists(edges_file)
)

############################################################
## 3) User-locked settings
############################################################
ms1_high_value <- "MS1_high"

# IL1R2 comparison within MS1_high (reference first)
il1r2_levels_within <- c("IL1R2_zero", "IL1R2_positive")  # positive vs zero

mother_keep <- c(
  "Innate Immune System",
  "Adaptive Immune System",
  "Cytokine Signaling in Immune system",
  "Hemostasis"
)

fdr_cutoff_union <- 0.05
bar_width <- 0.72

############################################################
## 4) Helper: load counts as integer matrix (FIXED)
############################################################
load_counts_int <- function(file) {
  cnt <- read.csv(file, check.names = FALSE)
  
  ens_raw   <- as.character(cnt[[1]])
  ens_clean <- sub("\\..*$", "", ens_raw)
  
  cnt[[1]] <- NULL
  rownames(cnt) <- ens_clean
  
  cnt <- as.matrix(cnt)
  
  storage.mode(cnt) <- "numeric"
  if (any(cnt < 0, na.rm = TRUE)) stop("Counts contain negative values.")
  if (any(is.na(cnt))) stop("Counts contain NA values.")
  
  cnt_int <- round(cnt)
  if (!all(cnt_int == cnt)) stop("Counts are not integer-like. DESeq2 requires raw integer counts.")
  storage.mode(cnt_int) <- "integer"
  
  cnt_int
}

############################################################
## 5) Helper: align phenotype to counts
############################################################
align_meta_counts <- function(cnt_int, studydata, group_col, levels_vec) {
  studydata <- studydata %>% dplyr::mutate(Patient_ID = as.character(Patient_ID))
  
  required_cols <- c("Patient_ID", group_col)
  missing_cols  <- setdiff(required_cols, colnames(studydata))
  if (length(missing_cols) > 0) {
    stop(paste0("Phenotype file missing columns: ", paste(missing_cols, collapse = ", ")))
  }
  
  missing_in_counts <- setdiff(studydata$Patient_ID, colnames(cnt_int))
  if (length(missing_in_counts) > 0) {
    stop(paste0(
      "These Patient_ID are in phenotype but not in counts: ",
      paste(head(missing_in_counts, 30), collapse = ", "),
      ifelse(length(missing_in_counts) > 30, " ...", "")
    ))
  }
  
  cnt_int <- cnt_int[, studydata$Patient_ID, drop = FALSE]
  stopifnot(all(colnames(cnt_int) == studydata$Patient_ID))
  
  meta <- studydata %>%
    dplyr::transmute(
      Patient_ID = Patient_ID,
      group = factor(.data[[group_col]], levels = levels_vec)
    )
  
  if (any(is.na(meta$group))) {
    stop(paste0(group_col, " has unexpected values. Expected levels: ", paste(levels_vec, collapse = ", ")))
  }
  
  list(cnt_int = cnt_int, meta = meta)
}

############################################################
## 6) Helper: get top1000 signature SYMBOLS (MS1, IL1R2)
############################################################
get_top1000_syms_ms1 <- function(sig_file, n = 1000) {
  sig <- read.delim(sig_file, stringsAsFactors = FALSE)
  req <- c("gennames", "MS1")
  if (!all(req %in% colnames(sig))) stop("MS1 signature file must contain columns: gennames, MS1")
  
  sig %>%
    dplyr::arrange(dplyr::desc(abs(MS1))) %>%
    dplyr::slice_head(n = n) %>%
    dplyr::pull(gennames) %>%
    as.character()
}

get_top1000_syms_il1r2 <- function(sig_file, n = 1000) {
  sig <- read.delim(sig_file, stringsAsFactors = FALSE, check.names = FALSE)
  
  gene_col  <- "NAME"
  score_col <- "IL1R2+_immature_neutrophils"
  
  if (!all(c(gene_col, score_col) %in% colnames(sig))) {
    stop(
      paste0(
        "IL1R2 signature file does not contain expected columns.\n",
        "Expected gene_col  = '", gene_col, "'\n",
        "Expected score_col = '", score_col, "'\n\n",
        "Found columns:\n",
        paste(colnames(sig), collapse = ", ")
      )
    )
  }
  
  sig %>%
    dplyr::transmute(
      gene  = as.character(.data[[gene_col]]),
      score = suppressWarnings(as.numeric(.data[[score_col]]))
    ) %>%
    dplyr::filter(!is.na(gene), gene != "", !is.na(score)) %>%
    dplyr::arrange(dplyr::desc(abs(score))) %>%
    dplyr::slice_head(n = n) %>%
    dplyr::pull(gene) %>%
    as.character()
}

############################################################
## 7) Helper: remove UNION signature genes (MS1 top1000 ∪ IL1R2 top1000)
##    RIGOROUS mapping: multiVals="list" + unlist (remove ALL mapped ENSEMBL)
############################################################
remove_union_signature_top1000 <- function(cnt_int, ms1_sig_file, il1r2_sig_file, n = 1000) {
  
  top_syms_ms1 <- get_top1000_syms_ms1(ms1_sig_file, n = n)
  top_syms_il1 <- get_top1000_syms_il1r2(il1r2_sig_file, n = n)
  
  top_syms_union <- unique(c(top_syms_ms1, top_syms_il1))
  
  top_ens_list <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys      = unique(top_syms_union),
    keytype   = "SYMBOL",
    column    = "ENSEMBL",
    multiVals = "list"
  )
  top_ens <- unique(na.omit(unlist(top_ens_list)))
  
  message("Union signature removal (ALL mappings): removing ", length(top_ens), " ENSEMBL genes (unique).")
  
  cnt_int[!rownames(cnt_int) %in% top_ens, , drop = FALSE]
}

############################################################
## 8) Helper: run DESeq2 + Reactome GSEA and return reactome_df
##    - break ties in stats + eps=0 (same as your 2nd script)
############################################################
run_deseq2_reactome <- function(cnt_int, meta, contrast_vec) {
  
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = cnt_int,
    colData   = meta %>% dplyr::select(group),
    design    = ~ group
  )
  
  keep <- rowSums(DESeq2::counts(dds) >= 10) >= 3
  dds  <- dds[keep, ]
  
  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::results(dds, contrast = contrast_vec)
  
  res_df <- as.data.frame(res)
  res_df$ENSEMBL <- rownames(res_df)
  res_df <- res_df %>% dplyr::filter(!is.na(stat))
  
  res_df$entrez <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys      = res_df$ENSEMBL,
    keytype   = "ENSEMBL",
    column    = "ENTREZID",
    multiVals = "first"
  )
  
  gene_df <- res_df %>%
    dplyr::filter(!is.na(entrez)) %>%
    dplyr::group_by(entrez) %>%
    dplyr::slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
  
  gene_list <- gene_df$stat
  names(gene_list) <- gene_df$entrez
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  # break ties deterministically (reduces fgsea warnings)
  if (any(duplicated(gene_list))) {
    set.seed(11)
    gene_list <- gene_list + runif(length(gene_list), min = -1e-7, max = 1e-7)
    gene_list <- sort(gene_list, decreasing = TRUE)
  }
  
  set.seed(11)
  gsea <- ReactomePA::gsePathway(
    geneList      = gene_list,
    organism      = "human",
    pvalueCutoff  = 1,
    pAdjustMethod = "BH",
    minGSSize     = 10,
    maxGSSize     = 10000000,
    verbose       = FALSE,
    nPermSimple   = 10000,
    eps           = 0
  )
  
  as.data.frame(gsea)
}

############################################################
## 9) Prepare Reactome backbone (Mother / Children)
############################################################
reactome_edges <- read.delim(edges_file) %>%
  dplyr::mutate(
    Mother   = stringr::str_trim(Mother),
    Children = stringr::str_trim(Children)
  )

if (!all(c("Mother","Children") %in% colnames(reactome_edges))) {
  stop("reactome_gl.txt must contain columns: Mother and Children")
}

reactome_edges <- reactome_edges %>%
  dplyr::filter(Mother %in% mother_keep)

mother_levels <- mother_keep[mother_keep %in% reactome_edges$Mother]

############################################################
## 10) Load phenotype once
############################################################
studydata <- readr::read_csv(pheno_file, show_col_types = FALSE) %>%
  dplyr::mutate(Patient_ID = as.character(Patient_ID))

############################################################
## 11) Subset: MS1_high only, then align IL1R2 groups to counts
############################################################
study_ms1_high <- studydata %>%
  dplyr::filter(MS1_group == ms1_high_value) %>%
  dplyr::filter(IL1R2_group %in% il1r2_levels_within) %>%
  dplyr::mutate(IL1R2_group = factor(IL1R2_group, levels = il1r2_levels_within))

message("MS1_high subset N = ", nrow(study_ms1_high))
message("IL1R2 within MS1_high counts:")
print(table(study_ms1_high$IL1R2_group))

cnt0 <- load_counts_int(counts_file)

aligned_il1_within <- align_meta_counts(
  cnt_int    = cnt0,
  studydata  = study_ms1_high,
  group_col  = "IL1R2_group",
  levels_vec = il1r2_levels_within
)

cnt_within  <- aligned_il1_within$cnt_int
meta_within <- aligned_il1_within$meta

############################################################
## 12) Circularity control: remove UNION top1000 (MS1 ∪ IL1R2)
############################################################
cnt_within_cc <- remove_union_signature_top1000(
  cnt_int        = cnt_within,
  ms1_sig_file   = ms1_sig_file,
  il1r2_sig_file = il1r2_sig_file,
  n              = 1000
)

############################################################
## 13) Run Reactome GSEA: IL1R2_positive vs IL1R2_zero within MS1_high
############################################################
reactome_df_within <- run_deseq2_reactome(
  cnt_int      = cnt_within_cc,
  meta         = meta_within,
  contrast_vec = c("group", "IL1R2_positive", "IL1R2_zero")
)

############################################################
## 14) Prepare GSEA table (global p.adjust + pvalue + NES)
############################################################
prep_gsea_tbl <- function(reactome_df, label) {
  reactome_df %>%
    dplyr::mutate(
      Description = stringr::str_trim(Description),
      NES         = as.numeric(NES),
      p.adjust    = as.numeric(p.adjust),
      pvalue      = as.numeric(pvalue)
    ) %>%
    dplyr::select(Description, NES, pvalue, p.adjust) %>%
    dplyr::mutate(analysis = label)
}

gsea_within <- prep_gsea_tbl(reactome_df_within, "IL1R2 within MS1_high")

############################################################
## 15) Child pathway selection: GLOBAL BH against backbone children
############################################################
all_children <- reactome_edges %>% dplyr::pull(Children) %>% unique()

master_children <- gsea_within %>%
  dplyr::filter(Description %in% all_children, !is.na(p.adjust), p.adjust < fdr_cutoff_union) %>%
  dplyr::pull(Description) %>%
  unique() %>%
  sort()

############################################################
## 16) Build backbone rows (Mother headers + selected Children)
############################################################
df_mothers <- tibble::tibble(
  parent_cat  = mother_levels,
  Description = mother_levels,
  is_mother   = TRUE
)

df_children <- reactome_edges %>%
  dplyr::filter(Children %in% master_children) %>%
  dplyr::transmute(
    parent_cat  = Mother,
    Description = Children,
    is_mother   = FALSE
  )

df_backbone <- purrr::map_dfr(mother_levels, function(m) {
  dplyr::bind_rows(
    df_mothers  %>% dplyr::filter(parent_cat == m),
    df_children %>% dplyr::filter(parent_cat == m)
  )
})

############################################################
## 17) Join GSEA results (GLOBAL BH) to backbone
############################################################
df_plot <- df_backbone %>%
  dplyr::left_join(
    gsea_within %>% dplyr::select(Description, NES, pvalue, p.adjust),
    by = "Description"
  ) %>%
  dplyr::mutate(
    analysis = "IL1R2 within MS1_high",
    
    y_label = ifelse(is_mother, paste0("**", parent_cat, "**"), Description),
    
    NES_plot = dplyr::case_when(
      is_mother & is.na(NES) ~ 0,
      TRUE ~ NES
    ),
    
    sig_global = dplyr::case_when(
      !is.na(p.adjust) & p.adjust < 0.05 ~ "FDR < 0.05",
      TRUE ~ "Not significant"
    )
  )

y_levels <- rev(unique(df_plot$y_label))
df_plot$y_label <- factor(df_plot$y_label, levels = y_levels)

nes_range <- range(df_plot$NES_plot, na.rm = TRUE)
nes_lim   <- max(abs(nes_range))
if (!is.finite(nes_lim) || nes_lim == 0) nes_lim <- 1

df_plot <- df_plot %>%
  dplyr::mutate(
    bar_fill = ifelse(!is.na(NES_plot) & sig_global == "FDR < 0.05", NES_plot, NA_real_)
  )

############################################################
## 18) Plot: single-panel BAR (GLOBAL BH) with desired theme
############################################################
p_within_globalBH <- ggplot2::ggplot(df_plot, ggplot2::aes(x = NES_plot, y = y_label)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 1.05, color = "black") +
  ggplot2::geom_col(
    ggplot2::aes(fill = bar_fill),
    width = bar_width,
    color = "black",
    linewidth = 0.25,
    na.rm = FALSE
  ) +
  ggplot2::scale_fill_gradient2(
    low  = "#4575b4",
    mid  = "white",
    high = "#d73027",
    limits = c(-nes_lim, nes_lim),
    name = "NES",
    na.value = "grey75"
  ) +
  ggplot2::theme_classic(base_size = BASE_SIZE, base_family = FONT_FAMILY) +
  ggplot2::theme(
    text         = ggplot2::element_text(family = FONT_FAMILY),
    axis.text.y  = ggtext::element_markdown(size = AXIS_TEXT, family = FONT_FAMILY),
    axis.text.x  = ggplot2::element_text(size = AXIS_TEXT, family = FONT_FAMILY),
    axis.title.x = ggplot2::element_text(size = AXIS_TITLE, family = FONT_FAMILY),
    axis.ticks.y = ggplot2::element_blank(),
    axis.line.y  = ggplot2::element_blank(),
    axis.line.x  = ggplot2::element_line(linewidth = LINE_AXIS, color = "black"),
    axis.ticks.x = ggplot2::element_line(linewidth = LINE_TICKS, color = "black"),
    panel.grid   = ggplot2::element_blank(),
    legend.position = "right",
    legend.title = ggplot2::element_text(size = AXIS_TITLE, family = FONT_FAMILY),
    legend.text  = ggplot2::element_text(size = AXIS_TEXT, family = FONT_FAMILY),
    axis.title.y = ggplot2::element_blank(),
    plot.margin  = ggplot2::margin(
      t = MARGINS_MM[1], r = MARGINS_MM[2], b = MARGINS_MM[3], l = MARGINS_MM[4],
      unit = "mm"
    )
  ) +
  labs(x = "Normalized Enrichment Score (NES)", y = NULL)

print(p_within_globalBH)

############################################################
## 19) Save
############################################################
dir.create("Figure", showWarnings = FALSE, recursive = TRUE)

ggsave(
  filename = "Figure/Reactome_IL1R2_bar_globalBH_within_MS1_high_unionSigRemoved.svg",
  plot     = p_within_globalBH,
  device   = svglite::svglite,
  width    = 12,
  height   = 8,
  units    = "in"
)
ggsave(
  filename    = "Figure/Reactome_IL1R2_bar_globalBH_within_MS1_high_unionSigRemoved.tiff",
  plot        = p_within_globalBH,
  device      = ragg::agg_tiff,
  width       = 12,
  height      = 8,
  units       = "in",
  dpi         = 600,
  compression = "lzw"
)
############################################################
## END
############################################################










#not exclude the signature genes####
############################################################
## Reactome GSEA: MS1 vs IL1R2+ SIDE-BY-SIDE BAR PLOT
## (NO circularity control: DO NOT remove signature genes)
## Everything else kept the same as your "signature-removed" script.
############################################################

############################################################
## 0) Clean start
############################################################
rm(list = ls())

############################################################
## 1) Packages
############################################################
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

pkgs_cran <- c(
  "dplyr","readr","tibble","stringr","tidyr","purrr",
  "ggplot2","scales","ggtext","patchwork"
)
pkgs_bioc <- c(
  "DESeq2","ReactomePA","org.Hs.eg.db","AnnotationDbi"
)

pacman::p_load(char = c(pkgs_cran, pkgs_bioc))

############################################################
## Global style (LOCKED)
############################################################
FONT_FAMILY <- "Arial"  # Font family used for all text
BASE_SIZE   <- 14       # Base font size (pt)
TITLE_SIZE  <- 14       # Plot title size (pt)
TAG_SIZE    <- 16       # Panel tag size (pt) - NOT bold
AXIS_TITLE  <- 13       # Axis title size (pt)
AXIS_TEXT   <- 11       # Axis tick label size (pt)

LINE_AXIS   <- 0.6      # Axis line width
LINE_TICKS  <- 0.6      # Tick line width
BIN_EDGE_W  <- 0.2      # Histogram bin edge width
MARGINS_MM  <- c(5.5, 5.5, 5.5, 5.5)  # Plot margins (t, r, b, l) in mm

############################################################
## 2) Inputs (files)
############################################################
counts_file    <- "Original_Data/Combined_counts_Amsterdam_my_cohort_allpatients.csv"
pheno_file     <- "Original_Data/studydata_with_ms1_il1r2.csv"
ms1_sig_file   <- "Original_Data/MS1_septiStates_150_Reyes_2020.txt"
il1r2_sig_file <- "Original_data/CIBERSORTx_Job16_ref_matrix_equal_inferred_phenoclasses.CIBERSORTx_Job16_ref_matrix_equal_inferred_refsample.bm.K999.txt"
edges_file     <- "Original_Data/reactome_gl.txt"

stopifnot(
  file.exists(counts_file),
  file.exists(pheno_file),
  file.exists(ms1_sig_file),
  file.exists(il1r2_sig_file),
  file.exists(edges_file)
)

############################################################
## 3) User-locked settings
############################################################
# Comparisons (locked)
ms1_levels   <- c("MS1_low", "MS1_high")                 # MS1_high vs MS1_low
il1r2_levels <- c("IL1R2_zero", "IL1R2_positive")        # IL1R2_positive vs IL1R2_zero

# Mother pathways (biological definition)
mother_keep <- c(
  "Innate Immune System",
  "Adaptive Immune System",
  "Cytokine Signaling in Immune system",
  "Hemostasis"
)

# GLOBAL BH threshold
fdr_cutoff_union <- 0.05

# Optional fallback (only when a mother would be empty in DISPLAY)
# (kept, but now "significance" still comes from global p.adjust; fallback is only for display)
fallback_n_max     <- 2
fallback_nominal_p <- 0.05

# Plot settings
bar_width <- 0.72

############################################################
## 4) Helper: load counts as integer matrix (FIXED)
############################################################
load_counts_int <- function(file) {
  cnt <- read.csv(file, check.names = FALSE)
  
  # First column = Ensembl IDs (possibly with version)
  ens_raw   <- as.character(cnt[[1]])
  ens_clean <- sub("\\..*$", "", ens_raw)
  
  cnt[[1]] <- NULL
  rownames(cnt) <- ens_clean
  
  # Force matrix
  cnt <- as.matrix(cnt)
  
  # Checks
  storage.mode(cnt) <- "numeric"
  if (any(cnt < 0, na.rm = TRUE)) stop("Counts contain negative values.")
  if (any(is.na(cnt))) stop("Counts contain NA values.")
  
  # DESeq2 requires integer counts
  cnt_int <- round(cnt)
  if (!all(cnt_int == cnt)) stop("Counts are not integer-like. DESeq2 requires raw integer counts.")
  storage.mode(cnt_int) <- "integer"
  
  cnt_int
}

############################################################
## 5) Helper: align phenotype to counts
############################################################
align_meta_counts <- function(cnt_int, studydata, group_col, levels_vec) {
  studydata <- studydata %>% dplyr::mutate(Patient_ID = as.character(Patient_ID))
  
  required_cols <- c("Patient_ID", group_col)
  missing_cols  <- setdiff(required_cols, colnames(studydata))
  if (length(missing_cols) > 0) {
    stop(paste0("Phenotype file missing columns: ", paste(missing_cols, collapse = ", ")))
  }
  
  missing_in_counts <- setdiff(studydata$Patient_ID, colnames(cnt_int))
  if (length(missing_in_counts) > 0) {
    stop(paste0(
      "These Patient_ID are in phenotype but not in counts: ",
      paste(head(missing_in_counts, 30), collapse = ", "),
      ifelse(length(missing_in_counts) > 30, " ...", "")
    ))
  }
  
  cnt_int <- cnt_int[, studydata$Patient_ID, drop = FALSE]
  stopifnot(all(colnames(cnt_int) == studydata$Patient_ID))
  
  meta <- studydata %>%
    dplyr::transmute(
      Patient_ID = Patient_ID,
      group = factor(.data[[group_col]], levels = levels_vec)
    )
  
  if (any(is.na(meta$group))) {
    stop(paste0(group_col, " has unexpected values. Expected levels: ", paste(levels_vec, collapse = ", ")))
  }
  
  list(cnt_int = cnt_int, meta = meta)
}

############################################################
## 6) Helper: run DESeq2 + Reactome GSEA and return reactome_df
##    Targeted update: break ties in stats + eps=0
############################################################
run_deseq2_reactome <- function(cnt_int, meta, contrast_vec) {
  
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = cnt_int,
    colData   = meta %>% dplyr::select(group),
    design    = ~ group
  )
  
  # Low count filtering
  keep <- rowSums(DESeq2::counts(dds) >= 10) >= 3
  dds  <- dds[keep, ]
  
  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::results(dds, contrast = contrast_vec)
  
  res_df <- as.data.frame(res)
  res_df$ENSEMBL <- rownames(res_df)
  res_df <- res_df %>% dplyr::filter(!is.na(stat))
  
  res_df$entrez <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys      = res_df$ENSEMBL,
    keytype   = "ENSEMBL",
    column    = "ENTREZID",
    multiVals = "first"
  )
  
  gene_df <- res_df %>%
    dplyr::filter(!is.na(entrez)) %>%
    dplyr::group_by(entrez) %>%
    dplyr::slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
  
  gene_list <- gene_df$stat
  names(gene_list) <- gene_df$entrez
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  # break ties deterministically
  if (any(duplicated(gene_list))) {
    set.seed(11)
    gene_list <- gene_list + runif(length(gene_list), min = -1e-7, max = 1e-7)
    gene_list <- sort(gene_list, decreasing = TRUE)
  }
  
  set.seed(11)
  gsea <- ReactomePA::gsePathway(
    geneList      = gene_list,
    organism      = "human",
    pvalueCutoff  = 1,
    pAdjustMethod = "BH",
    minGSSize     = 10,
    maxGSSize     = 10000000,
    verbose       = FALSE,
    nPermSimple   = 15000,
    eps           = 0
  )
  
  as.data.frame(gsea)
}

############################################################
## 7) Prepare Reactome backbone (Mother / Children)
############################################################
reactome_edges <- read.delim(edges_file) %>%
  dplyr::mutate(
    Mother   = stringr::str_trim(Mother),
    Children = stringr::str_trim(Children)
  )

if (!all(c("Mother","Children") %in% colnames(reactome_edges))) {
  stop("reactome_gl.txt must contain columns: Mother and Children")
}

reactome_edges <- reactome_edges %>%
  dplyr::filter(Mother %in% mother_keep)

# LOCK mother order exactly as mother_keep
mother_levels <- mother_keep[mother_keep %in% reactome_edges$Mother]

############################################################
## 8) Load phenotype once
############################################################
studydata <- readr::read_csv(pheno_file, show_col_types = FALSE) %>%
  dplyr::mutate(Patient_ID = as.character(Patient_ID))

############################################################
## 9) Run MS1 Reactome (NO signature removal)
############################################################
cnt0 <- load_counts_int(counts_file)

aligned_ms1 <- align_meta_counts(cnt0, studydata, group_col = "MS1_group", levels_vec = ms1_levels)
cnt_ms1     <- aligned_ms1$cnt_int
meta_ms1    <- aligned_ms1$meta

# NO circularity control here
reactome_df_ms1 <- run_deseq2_reactome(
  cnt_int      = cnt_ms1,
  meta         = meta_ms1,
  contrast_vec = c("group", "MS1_high", "MS1_low")
)

############################################################
## 10) Run IL1R2 Reactome (NO signature removal)
############################################################
aligned_il1 <- align_meta_counts(cnt0, studydata, group_col = "IL1R2_group", levels_vec = il1r2_levels)
cnt_il1     <- aligned_il1$cnt_int
meta_il1    <- aligned_il1$meta

# NO circularity control here
reactome_df_il1r2 <- run_deseq2_reactome(
  cnt_int      = cnt_il1,
  meta         = meta_il1,
  contrast_vec = c("group", "IL1R2_positive", "IL1R2_zero")
)

############################################################
## 11) Prepare GSEA tables (global p.adjust + pvalue + NES)
############################################################
prep_gsea_tbl <- function(reactome_df, label) {
  reactome_df %>%
    dplyr::mutate(
      Description = stringr::str_trim(Description),
      NES         = as.numeric(NES),
      p.adjust    = as.numeric(p.adjust),
      pvalue      = as.numeric(pvalue)
    ) %>%
    dplyr::select(Description, NES, pvalue, p.adjust) %>%
    dplyr::mutate(analysis = label)
}

gsea_ms1 <- prep_gsea_tbl(reactome_df_ms1,   "MS1")
gsea_il1 <- prep_gsea_tbl(reactome_df_il1r2, "IL1R2")

############################################################
## 12) Child pathway selection: UNION rule using GLOBAL BH (p.adjust)
############################################################
all_children <- reactome_edges %>% dplyr::pull(Children) %>% unique()

sig_children_ms1 <- gsea_ms1 %>%
  dplyr::filter(Description %in% all_children, !is.na(p.adjust), p.adjust < fdr_cutoff_union) %>%
  dplyr::pull(Description)

sig_children_il1 <- gsea_il1 %>%
  dplyr::filter(Description %in% all_children, !is.na(p.adjust), p.adjust < fdr_cutoff_union) %>%
  dplyr::pull(Description)

master_children <- sort(unique(c(sig_children_ms1, sig_children_il1)))

############################################################
## 13) Optional fallback: if a Mother would be empty, add up to 2 children
############################################################
fallback_children <- c()

for (m in mother_levels) {
  children_m <- reactome_edges %>% dplyr::filter(Mother == m) %>% dplyr::pull(Children) %>% unique()
  shown_m    <- intersect(children_m, master_children)
  
  if (length(shown_m) == 0) {
    
    cand_ms1 <- gsea_ms1 %>%
      dplyr::filter(Description %in% children_m, !is.na(pvalue), pvalue < fallback_nominal_p) %>%
      dplyr::arrange(pvalue) %>%
      dplyr::pull(Description)
    
    cand_il1 <- gsea_il1 %>%
      dplyr::filter(Description %in% children_m, !is.na(pvalue), pvalue < fallback_nominal_p) %>%
      dplyr::arrange(pvalue) %>%
      dplyr::pull(Description)
    
    cand_union <- unique(c(cand_ms1, cand_il1))
    
    if (length(cand_union) > 0) {
      fallback_children <- c(fallback_children, head(cand_union, fallback_n_max))
    }
  }
}

fallback_children <- sort(unique(fallback_children))
master_children2  <- sort(unique(c(master_children, fallback_children)))

############################################################
## 14) Build backbone rows (Mother headers + selected Children)
############################################################
df_mothers <- tibble::tibble(
  parent_cat  = mother_levels,
  Description = mother_levels,
  is_mother   = TRUE
)

df_children <- reactome_edges %>%
  dplyr::filter(Children %in% master_children2) %>%
  dplyr::transmute(
    parent_cat  = Mother,
    Description = Children,
    is_mother   = FALSE
  )

df_backbone <- purrr::map_dfr(mother_levels, function(m) {
  dplyr::bind_rows(
    df_mothers  %>% dplyr::filter(parent_cat == m),
    df_children %>% dplyr::filter(parent_cat == m)
  )
})

############################################################
## 15) Join GSEA results (GLOBAL BH) to backbone
############################################################
join_global <- function(df_backbone, gsea_tbl_one, analysis_label) {
  
  out <- df_backbone %>%
    dplyr::left_join(
      gsea_tbl_one %>% dplyr::select(Description, NES, pvalue, p.adjust),
      by = "Description"
    ) %>%
    dplyr::mutate(
      analysis = analysis_label,
      y_label = ifelse(is_mother, paste0("**", parent_cat, "**"), Description),
      NES_plot = dplyr::case_when(
        is_mother & is.na(NES) ~ 0,
        TRUE ~ NES
      ),
      sig_global = dplyr::case_when(
        !is.na(p.adjust) & p.adjust < 0.05 ~ "FDR < 0.05",
        TRUE ~ "Not significant"
      )
    )
  
  out
}

df_ms1 <- join_global(df_backbone, gsea_ms1, "MS1")
df_il1 <- join_global(df_backbone, gsea_il1, "IL1R2")

df_both <- dplyr::bind_rows(df_ms1, df_il1)
df_both$analysis <- factor(df_both$analysis, levels = c("MS1", "IL1R2"))

y_levels <- rev(unique(df_both$y_label))
df_both$y_label <- factor(df_both$y_label, levels = y_levels)

nes_range <- range(df_both$NES_plot, na.rm = TRUE)
nes_lim   <- max(abs(nes_range))
if (!is.finite(nes_lim) || nes_lim == 0) nes_lim <- 1

df_both <- df_both %>%
  dplyr::mutate(
    bar_fill = ifelse(!is.na(NES_plot) & sig_global == "FDR < 0.05", NES_plot, NA_real_)
  )

############################################################
## 16) Plot: side-by-side BAR (GLOBAL BH) with LOCKED style
############################################################
p_main_globalBH <- ggplot2::ggplot(df_both, ggplot2::aes(x = NES_plot, y = y_label)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 1.05, color = "black") +
  ggplot2::geom_col(
    ggplot2::aes(fill = bar_fill),
    width = bar_width,
    color = "black",
    linewidth = 0.25,
    na.rm = FALSE
  ) +
  ggplot2::facet_wrap(~analysis, nrow = 1, scales = "free_x") +
  ggplot2::scale_fill_gradient2(
    low  = "#4575b4",
    mid  = "white",
    high = "#d73027",
    limits = c(-nes_lim, nes_lim),
    name = "NES",
    na.value = "grey75"
  ) +
  ggplot2::theme_classic(base_size = BASE_SIZE, base_family = FONT_FAMILY) +
  ggplot2::theme(
    text               = ggplot2::element_text(family = FONT_FAMILY),
    plot.title         = ggplot2::element_text(size = TITLE_SIZE, family = FONT_FAMILY),
    axis.text.y        = ggtext::element_markdown(size = AXIS_TEXT, family = FONT_FAMILY),
    axis.text.x        = ggplot2::element_text(size = AXIS_TEXT, family = FONT_FAMILY),
    axis.title.x       = ggplot2::element_text(size = AXIS_TITLE, family = FONT_FAMILY),
    axis.title.y       = ggplot2::element_text(size = AXIS_TITLE, family = FONT_FAMILY),
    axis.ticks.y       = ggplot2::element_blank(),
    axis.ticks.x       = ggplot2::element_line(linewidth = LINE_TICKS, color = "black"),
    axis.line.y        = ggplot2::element_blank(),
    axis.line.x        = ggplot2::element_line(linewidth = LINE_AXIS, color = "black"),
    panel.grid         = ggplot2::element_blank(),
    strip.background   = ggplot2::element_blank(),
    strip.text         = ggplot2::element_text(face = "bold", size = BASE_SIZE, family = FONT_FAMILY),
    legend.position    = "right",
    legend.title       = ggplot2::element_text(size = AXIS_TITLE, family = FONT_FAMILY),
    legend.text        = ggplot2::element_text(size = AXIS_TEXT, family = FONT_FAMILY),
    axis.title         = ggplot2::element_blank(),
    plot.margin        = ggplot2::margin(
      t = MARGINS_MM[1], r = MARGINS_MM[2], b = MARGINS_MM[3], l = MARGINS_MM[4],
      unit = "mm"
    )
  ) +
  labs(x = "Normalized Enrichment Score (NES)", y = NULL)

print(p_main_globalBH)

############################################################
## 17) OPTIONAL: Save
############################################################
if (!requireNamespace("svglite", quietly = TRUE)) install.packages("svglite")
library(svglite)

dir.create("Figure", showWarnings = FALSE, recursive = TRUE)

ggsave(
  filename = "Figure/Reactome_MS1_vs_IL1R2_bar_globalBH_main_NOsignatureRemoval.svg",
  plot     = p_main_globalBH,
  device   = svglite::svglite,
  width    = 12,
  height   = 8,
  units    = "in"
)

# not exclude the signature genes within MS1 high group######
#NO exclusion of MS1/IL1R2 signatures within MS1_high####
############################################################
## 0) Clean start
############################################################
rm(list = ls())

############################################################
## 1) Packages
############################################################
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

pkgs_cran <- c(
  "dplyr","readr","tibble","stringr","tidyr","purrr",
  "ggplot2","scales","ggtext"
)
pkgs_bioc <- c(
  "DESeq2","ReactomePA","org.Hs.eg.db","AnnotationDbi"
)

pacman::p_load(char = c(pkgs_cran, pkgs_bioc))

############################################################
## Global style (LOCKED)
############################################################
FONT_FAMILY <- "Arial"  # Font family used for all text
BASE_SIZE   <- 14       # Base font size (pt)
TITLE_SIZE  <- 14       # Plot title size (pt)
TAG_SIZE    <- 16       # Panel tag size (pt) - NOT bold
AXIS_TITLE  <- 13       # Axis title size (pt)
AXIS_TEXT   <- 11       # Axis tick label size (pt)

LINE_AXIS   <- 0.6      # Axis line width
LINE_TICKS  <- 0.6      # Tick line width
BIN_EDGE_W  <- 0.2      # Histogram bin edge width
MARGINS_MM  <- c(5.5, 5.5, 5.5, 5.5)  # Plot margins (t, r, b, l) in mm

############################################################
## 2) Inputs (files)
############################################################
counts_file    <- "Original_Data/Combined_counts_Amsterdam_my_cohort_allpatients.csv"
pheno_file     <- "Original_Data/studydata_with_ms1_il1r2.csv"
ms1_sig_file   <- "Original_Data/MS1_septiStates_150_Reyes_2020.txt"
il1r2_sig_file <- "Original_data/CIBERSORTx_Job16_ref_matrix_equal_inferred_phenoclasses.CIBERSORTx_Job16_ref_matrix_equal_inferred_refsample.bm.K999.txt"
edges_file     <- "Original_Data/reactome_gl.txt"

stopifnot(
  file.exists(counts_file),
  file.exists(pheno_file),
  file.exists(ms1_sig_file),
  file.exists(il1r2_sig_file),
  file.exists(edges_file)
)

############################################################
## 3) User-locked settings
############################################################
# MS1 high subset
ms1_high_value <- "MS1_high"

# IL1R2 comparison within MS1_high (reference first)
il1r2_levels_within <- c("IL1R2_zero", "IL1R2_positive")  # IL1R2_positive vs IL1R2_zero

# Mother pathways (biological definition) - LOCKED order
mother_keep <- c(
  "Innate Immune System",
  "Adaptive Immune System",
  "Cytokine Signaling in Immune system",
  "Hemostasis"
)

# GLOBAL BH threshold for child display
fdr_cutoff_union <- 0.05

# Plot settings
bar_width <- 0.72

############################################################
## 4) Helper: load counts as integer matrix (FIXED)
############################################################
load_counts_int <- function(file) {
  cnt <- read.csv(file, check.names = FALSE)
  
  # First column = Ensembl IDs (possibly with version)
  ens_raw   <- as.character(cnt[[1]])
  ens_clean <- sub("\\..*$", "", ens_raw)
  
  cnt[[1]] <- NULL
  rownames(cnt) <- ens_clean
  
  # Force matrix
  cnt <- as.matrix(cnt)
  
  # Checks
  storage.mode(cnt) <- "numeric"
  if (any(cnt < 0, na.rm = TRUE)) stop("Counts contain negative values.")
  if (any(is.na(cnt))) stop("Counts contain NA values.")
  
  # DESeq2 requires integer counts
  cnt_int <- round(cnt)
  if (!all(cnt_int == cnt)) stop("Counts are not integer-like. DESeq2 requires raw integer counts.")
  storage.mode(cnt_int) <- "integer"
  
  cnt_int
}

############################################################
## 5) Helper: align phenotype to counts
############################################################
align_meta_counts <- function(cnt_int, studydata, group_col, levels_vec) {
  studydata <- studydata %>% dplyr::mutate(Patient_ID = as.character(Patient_ID))
  
  required_cols <- c("Patient_ID", group_col)
  missing_cols  <- setdiff(required_cols, colnames(studydata))
  if (length(missing_cols) > 0) {
    stop(paste0("Phenotype file missing columns: ", paste(missing_cols, collapse = ", ")))
  }
  
  missing_in_counts <- setdiff(studydata$Patient_ID, colnames(cnt_int))
  if (length(missing_in_counts) > 0) {
    stop(paste0(
      "These Patient_ID are in phenotype but not in counts: ",
      paste(head(missing_in_counts, 30), collapse = ", "),
      ifelse(length(missing_in_counts) > 30, " ...", "")
    ))
  }
  
  cnt_int <- cnt_int[, studydata$Patient_ID, drop = FALSE]
  stopifnot(all(colnames(cnt_int) == studydata$Patient_ID))
  
  meta <- studydata %>%
    dplyr::transmute(
      Patient_ID = Patient_ID,
      group = factor(.data[[group_col]], levels = levels_vec)
    )
  
  if (any(is.na(meta$group))) {
    stop(paste0(group_col, " has unexpected values. Expected levels: ", paste(levels_vec, collapse = ", ")))
  }
  
  list(cnt_int = cnt_int, meta = meta)
}

############################################################
## 6) Helper: run DESeq2 + Reactome GSEA and return reactome_df
## (kept identical to your current version)
############################################################
run_deseq2_reactome <- function(cnt_int, meta, contrast_vec) {
  
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = cnt_int,
    colData   = meta %>% dplyr::select(group),
    design    = ~ group
  )
  
  # Low count filtering
  keep <- rowSums(DESeq2::counts(dds) >= 10) >= 3
  dds  <- dds[keep, ]
  
  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::results(dds, contrast = contrast_vec)
  
  res_df <- as.data.frame(res)
  res_df$ENSEMBL <- rownames(res_df)
  res_df <- res_df %>% dplyr::filter(!is.na(stat))
  
  res_df$entrez <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys      = res_df$ENSEMBL,
    keytype   = "ENSEMBL",
    column    = "ENTREZID",
    multiVals = "first"
  )
  
  gene_df <- res_df %>%
    dplyr::filter(!is.na(entrez)) %>%
    dplyr::group_by(entrez) %>%
    dplyr::slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
  
  gene_list <- gene_df$stat
  names(gene_list) <- gene_df$entrez
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  set.seed(11)
  gsea <- ReactomePA::gsePathway(
    geneList      = gene_list,
    organism      = "human",
    pvalueCutoff  = 1,
    pAdjustMethod = "BH",    # Global BH
    minGSSize     = 10,
    maxGSSize     = 10000000,
    verbose       = FALSE,
    nPermSimple   = 10000
  )
  
  as.data.frame(gsea)
}

############################################################
## 7) Prepare Reactome backbone (Mother / Children)
############################################################
reactome_edges <- read.delim(edges_file) %>%
  dplyr::mutate(
    Mother   = stringr::str_trim(Mother),
    Children = stringr::str_trim(Children)
  )

if (!all(c("Mother","Children") %in% colnames(reactome_edges))) {
  stop("reactome_gl.txt must contain columns: Mother and Children")
}

reactome_edges <- reactome_edges %>%
  dplyr::filter(Mother %in% mother_keep)

# LOCK mother order EXACTLY as mother_keep
mother_levels <- mother_keep[mother_keep %in% reactome_edges$Mother]

############################################################
## 8) Load phenotype once
############################################################
studydata <- readr::read_csv(pheno_file, show_col_types = FALSE) %>%
  dplyr::mutate(Patient_ID = as.character(Patient_ID))

############################################################
## 9) Subset: MS1_high only, then align IL1R2 groups to counts
############################################################
study_ms1_high <- studydata %>%
  dplyr::filter(MS1_group == ms1_high_value) %>%
  dplyr::filter(IL1R2_group %in% il1r2_levels_within) %>%
  dplyr::mutate(IL1R2_group = factor(IL1R2_group, levels = il1r2_levels_within))

message("MS1_high subset N = ", nrow(study_ms1_high))
message("IL1R2 within MS1_high counts:")
print(table(study_ms1_high$IL1R2_group))

cnt0 <- load_counts_int(counts_file)

aligned_il1_within <- align_meta_counts(
  cnt_int     = cnt0,
  studydata   = study_ms1_high,
  group_col   = "IL1R2_group",
  levels_vec  = il1r2_levels_within
)

cnt_within  <- aligned_il1_within$cnt_int
meta_within <- aligned_il1_within$meta

############################################################
## 10) (CHANGED) NO circularity control: DO NOT remove signature genes
############################################################
cnt_within_cc <- cnt_within

############################################################
## 11) Run Reactome GSEA: IL1R2_positive vs IL1R2_zero within MS1_high
############################################################
reactome_df_within <- run_deseq2_reactome(
  cnt_int      = cnt_within_cc,
  meta         = meta_within,
  contrast_vec = c("group", "IL1R2_positive", "IL1R2_zero")
)

############################################################
## 12) Prepare GSEA table (global p.adjust + pvalue + NES)
############################################################
prep_gsea_tbl <- function(reactome_df, label) {
  reactome_df %>%
    dplyr::mutate(
      Description = stringr::str_trim(Description),
      NES         = as.numeric(NES),
      p.adjust    = as.numeric(p.adjust),
      pvalue      = as.numeric(pvalue)
    ) %>%
    dplyr::select(Description, NES, pvalue, p.adjust) %>%
    dplyr::mutate(analysis = label)
}

gsea_within <- prep_gsea_tbl(reactome_df_within, "IL1R2 within MS1_high")

############################################################
## 13) Child pathway selection: GLOBAL BH against backbone children
############################################################
all_children <- reactome_edges %>% dplyr::pull(Children) %>% unique()

master_children <- gsea_within %>%
  dplyr::filter(Description %in% all_children, !is.na(p.adjust), p.adjust < fdr_cutoff_union) %>%
  dplyr::pull(Description) %>%
  unique() %>%
  sort()

############################################################
## 14) Build backbone rows (Mother headers + selected Children)
############################################################
df_mothers <- tibble::tibble(
  parent_cat  = mother_levels,
  Description = mother_levels,
  is_mother   = TRUE
)

df_children <- reactome_edges %>%
  dplyr::filter(Children %in% master_children) %>%
  dplyr::transmute(
    parent_cat  = Mother,
    Description = Children,
    is_mother   = FALSE
  )

df_backbone <- purrr::map_dfr(mother_levels, function(m) {
  dplyr::bind_rows(
    df_mothers  %>% dplyr::filter(parent_cat == m),
    df_children %>% dplyr::filter(parent_cat == m)
  )
})

############################################################
## 15) Join GSEA results (GLOBAL BH) to backbone
############################################################
df_plot <- df_backbone %>%
  dplyr::left_join(
    gsea_within %>% dplyr::select(Description, NES, pvalue, p.adjust),
    by = "Description"
  ) %>%
  dplyr::mutate(
    analysis = "IL1R2 within MS1_high",
    y_label  = ifelse(is_mother, paste0("**", parent_cat, "**"), Description),
    NES_plot = dplyr::case_when(
      is_mother & is.na(NES) ~ 0,
      TRUE ~ NES
    ),
    sig_global = dplyr::case_when(
      !is.na(p.adjust) & p.adjust < 0.05 ~ "FDR < 0.05",
      TRUE ~ "Not significant"
    )
  )

y_levels <- rev(unique(df_plot$y_label))
df_plot$y_label <- factor(df_plot$y_label, levels = y_levels)

nes_range <- range(df_plot$NES_plot, na.rm = TRUE)
nes_lim   <- max(abs(nes_range))
if (!is.finite(nes_lim) || nes_lim == 0) nes_lim <- 1

df_plot <- df_plot %>%
  dplyr::mutate(
    bar_fill = ifelse(!is.na(NES_plot) & sig_global == "FDR < 0.05", NES_plot, NA_real_)
  )

############################################################
## 16) Plot: single-panel BAR (GLOBAL BH) with desired theme
############################################################
p_within_globalBH <- ggplot2::ggplot(df_plot, ggplot2::aes(x = NES_plot, y = y_label)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 1.05, color = "black") +
  ggplot2::geom_col(
    ggplot2::aes(fill = bar_fill),
    width = bar_width,
    color = "black",
    linewidth = 0.25,
    na.rm = FALSE
  ) +
  ggplot2::scale_fill_gradient2(
    low  = "#4575b4",
    mid  = "white",
    high = "#d73027",
    limits = c(-nes_lim, nes_lim),
    name = "NES",
    na.value = "grey75"
  ) +
  ggplot2::theme_classic(base_size = BASE_SIZE, base_family = FONT_FAMILY) +
  ggplot2::theme(
    text         = ggplot2::element_text(family = FONT_FAMILY),
    axis.text.y  = ggtext::element_markdown(size = AXIS_TEXT, family = FONT_FAMILY),
    axis.text.x  = ggplot2::element_text(size = AXIS_TEXT, family = FONT_FAMILY),
    axis.title.x = ggplot2::element_text(size = AXIS_TITLE, family = FONT_FAMILY),
    axis.ticks.y = ggplot2::element_blank(),
    axis.line.y  = ggplot2::element_blank(),
    axis.line.x  = ggplot2::element_line(linewidth = LINE_AXIS, color = "black"),
    axis.ticks.x = ggplot2::element_line(linewidth = LINE_TICKS, color = "black"),
    panel.grid   = ggplot2::element_blank(),
    legend.position = "right",
    legend.title = ggplot2::element_text(size = AXIS_TITLE, family = FONT_FAMILY),
    legend.text  = ggplot2::element_text(size = AXIS_TEXT, family = FONT_FAMILY),
    axis.title.y = ggplot2::element_blank(),
    plot.margin  = ggplot2::margin(
      t = MARGINS_MM[1], r = MARGINS_MM[2], b = MARGINS_MM[3], l = MARGINS_MM[4],
      unit = "mm"
    )
  ) +
  labs(x = "Normalized Enrichment Score (NES)", y = NULL)

print(p_within_globalBH)

############################################################
## 17) OPTIONAL: Save (only if svglite is available)
############################################################
if (!requireNamespace("svglite", quietly = TRUE)) install.packages("svglite")
library(svglite)

dir.create("Figure", showWarnings = FALSE, recursive = TRUE)

ggsave(
  filename = "Figure/Reactome_IL1R2_within_MS1_high_globalBH_NOsignatureRemoval.svg",
  plot     = p_within_globalBH,
  device   = svglite::svglite,
  width    = 12,
  height   = 8,
  units    = "in"
)
ggsave(
  filename    = "Figure/Reactome_IL1R2_within_MS1_high_globalBH_NOsignatureRemoval.tiff",
  plot        = p_within_globalBH,
  device      = ragg::agg_tiff,
  width       = 12,
  height      = 8,
  units       = "in",
  dpi         = 600,
  compression = "lzw"
)


