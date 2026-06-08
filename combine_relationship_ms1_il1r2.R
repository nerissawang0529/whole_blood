############################################################
## Clean workspace
############################################################
rm(list = ls())

############################################################
## Packages
############################################################
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(rms)
  library(patchwork)
  library(svglite)
  library(mgcv)
  library(readr)
  library(stringr)
  library(forcats)
  library(rstatix)
  library(purrr)
  library(grid)
  library(ggtext)
})

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
## Global style (LOCKED / unified)
############################################################
FONT_FAMILY <- "Arial"   # if Arial is not available, use "sans"
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
## Colors (LOCKED)
############################################################
col_ms1   <- "#2F4B7C"
col_il1r2 <- "#A23B3B"

############################################################
## Visual intensity
############################################################
ALPHA_POINT_MS1   <- 0.75
ALPHA_FILL_MS1    <- 0.22
ALPHA_RIBBON_MS1  <- 0.05

ALPHA_POINT_IL1R2 <- 0.75
ALPHA_BOX_FILL    <- 0.28

POINT_SIZE <- 1.9
LINE_MAIN  <- 1.25
STAR_SIZE  <- 4.2

############################################################
## Unified themes
############################################################
base_theme <- theme_classic(base_size = BASE_SIZE, base_family = FONT_FAMILY) +
  theme(
    text        = element_text(family = FONT_FAMILY),
    plot.title  = element_text(size = TITLE_SIZE, face = "plain", hjust = 0),
    plot.tag    = element_text(size = TAG_SIZE, face = "plain", family = FONT_FAMILY),
    axis.title  = element_text(size = AXIS_TITLE, family = FONT_FAMILY),
    axis.text   = element_text(size = AXIS_TEXT, family = FONT_FAMILY),
    axis.line   = element_line(linewidth = LINE_AXIS),
    axis.ticks  = element_line(linewidth = LINE_TICKS),
    plot.margin = margin(MARGINS_MM[1], MARGINS_MM[2], MARGINS_MM[3], MARGINS_MM[4], unit = "mm")
  )

forest_theme <- base_theme +
  theme(
    axis.line.y = element_blank(),
    legend.position = "none",
    axis.text.y = element_text(size = AXIS_TEXT, family = FONT_FAMILY),
    axis.text.x = element_text(size = AXIS_TEXT, family = FONT_FAMILY),
    axis.title.x = element_text(size = AXIS_TITLE, family = FONT_FAMILY),
    plot.margin = margin(5.5, 12, 6, 5.5, unit = "mm")
  )

reactome_theme <- base_theme +
  theme(
    axis.text.y  = ggtext::element_markdown(size = AXIS_TEXT - 0.5, family = FONT_FAMILY),
    axis.text.x  = element_text(size = AXIS_TEXT, family = FONT_FAMILY),
    axis.title.x = element_text(size = AXIS_TITLE, family = FONT_FAMILY),
    axis.ticks.y = element_blank(),
    axis.line.y  = element_blank(),
    panel.grid   = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = AXIS_TITLE, family = FONT_FAMILY),
    legend.text  = element_text(size = AXIS_TEXT, family = FONT_FAMILY),
    axis.title.y = element_blank(),
    plot.margin  = margin(5.5, 5.5, 6, 1.5, unit = "mm")
  )

############################################################
## Load data
############################################################
studydata <- read.csv("Original_data/studydata_with_ms1_il1r2.csv")
wb_long   <- read_csv("Original_data/whole_blood_stimuli_long.csv", show_col_types = FALSE)

############################################################
## Prepare core variables for panels A and B
############################################################
df <- studydata %>%
  transmute(
    MS1       = as.numeric(MS1_percent),
    IL1R2     = as.numeric(IL1R2_percent),
    IL1R2_pos = ifelse(as.numeric(IL1R2_percent) > 0, 1, 0),
    PSI       = suppressWarnings(as.numeric(PSI_new)),
    MEWS      = suppressWarnings(as.numeric(MEWS_score))
  ) %>%
  filter(!is.na(MS1), !is.na(IL1R2))

############################################################
## Build IL1R2 strata: 0 + sextiles among positives
############################################################
df <- df %>% mutate(IL1R2_group = NA_character_)
df$IL1R2_group[df$IL1R2 == 0] <- "0"

pos_idx <- which(df$IL1R2 > 0)
df$IL1R2_group[pos_idx] <- as.character(ntile(df$IL1R2[pos_idx], 6))
df$IL1R2_group <- factor(df$IL1R2_group, levels = c("0","1","2","3","4","5","6"))

df <- df %>% mutate(IL1R2_group_num = as.numeric(as.character(IL1R2_group)))

############################################################
## Helper: Spearman plot (MS1 on x) + GAM smooth
############################################################
spearman_plot <- function(data, xvar, yvar, xlab, ylab, point_col, integer_breaks = NULL) {
  sp <- cor.test(data[[xvar]], data[[yvar]], method = "spearman", use = "complete.obs")
  lab <- paste0("Spearman \u03C1 = ", round(sp$estimate, 2), ", p = ", signif(sp$p.value, 2))
  
  p <- ggplot(data, aes(.data[[xvar]], .data[[yvar]])) +
    geom_point(color = point_col, alpha = ALPHA_POINT_MS1, size = POINT_SIZE) +
    geom_smooth(
      method  = "gam",
      formula = y ~ s(x, bs = "cs"),
      se      = TRUE,
      color   = "black",
      fill    = point_col,
      alpha   = ALPHA_FILL_MS1,
      linewidth = 1.1
    ) +
    annotate(
      "text",
      x = Inf, y = Inf, label = lab,
      hjust = 1.1, vjust = 1.4, size = 3.8,
      family = FONT_FAMILY
    ) +
    labs(x = xlab, y = ylab) +
    base_theme
  
  if (!is.null(integer_breaks)) {
    p <- p + scale_y_continuous(breaks = integer_breaks)
  }
  p
}

############################################################
## Helper: Box + jitter + KW + GAM curve ONLY for 1..6
############################################################
box_jitter_kw_with_curve <- function(data, yvar, xlab, ylab, kw_p, integer_breaks = NULL) {
  
  df2 <- data %>%
    filter(!is.na(.data[[yvar]]), !is.na(IL1R2_group_num))
  
  df_for_curve <- df2 %>% filter(IL1R2_group_num >= 1)
  
  gam_formula <- as.formula(paste0(yvar, " ~ s(IL1R2_group_num, k = 4)"))
  gam_fit <- mgcv::gam(gam_formula, data = df_for_curve)
  
  new_curve <- data.frame(IL1R2_group_num = seq(1, 6, length.out = 200))
  pred <- predict(gam_fit, newdata = new_curve, se.fit = TRUE)
  
  curve_df <- new_curve %>%
    mutate(
      fit  = pred$fit,
      low  = pred$fit - 1.96 * pred$se.fit,
      high = pred$fit + 1.96 * pred$se.fit
    )
  
  p <- ggplot(df2, aes(x = IL1R2_group_num, y = .data[[yvar]])) +
    geom_boxplot(
      aes(group = IL1R2_group_num),
      width = 0.70,
      fill = col_il1r2,
      alpha = ALPHA_BOX_FILL,
      color = "black",
      linewidth = 0.6,
      outlier.shape = NA
    ) +
    geom_jitter(
      width = 0.18,
      size = POINT_SIZE,
      alpha = ALPHA_POINT_IL1R2,
      color = col_il1r2
    ) +
    geom_ribbon(
      data = curve_df,
      inherit.aes = FALSE,
      aes(x = IL1R2_group_num, ymin = low, ymax = high),
      fill = col_il1r2,
      alpha = 0.20
    ) +
    geom_line(
      data = curve_df,
      inherit.aes = FALSE,
      aes(x = IL1R2_group_num, y = fit),
      color = col_il1r2,
      linewidth = 1.1
    ) +
    scale_x_continuous(
      breaks = 0:6,
      labels = 0:6,
      limits = c(-0.5, 6.5)
    ) +
    annotate(
      "text",
      x = Inf, y = Inf,
      label = paste0("Kruskal–Wallis p = ", signif(kw_p, 2)),
      hjust = 1.1, vjust = 1.4,
      size = 3.8,
      family = FONT_FAMILY
    ) +
    labs(x = xlab, y = ylab) +
    base_theme
  
  if (!is.null(integer_breaks)) {
    p <- p + scale_y_continuous(breaks = integer_breaks)
  }
  p
}

############################################################
## Figure A and B
############################################################
df$IL1R2_pos <- factor(df$IL1R2_pos)

dd <- datadist(df); options(datadist = "dd")
fit <- lrm(IL1R2_pos ~ rcs(MS1, 3), data = df)
pred <- Predict(fit, MS1, fun = plogis)

p_A1 <- ggplot(pred, aes(MS1, yhat)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              fill = col_ms1, alpha = ALPHA_RIBBON_MS1) +
  geom_line(color = col_ms1, linewidth = LINE_MAIN) +
  labs(
    x = "Percentage of MS1 cells",
    y = "Predicted probability of IL1R2+ positivity"
  ) +
  base_theme

kw_ms1 <- kruskal.test(MS1 ~ IL1R2_group, data = df)

p_B1 <- box_jitter_kw_with_curve(
  data = df,
  yvar = "MS1",
  xlab = "IL1R2 strata",
  ylab = "Percentage of MS1 cells",
  kw_p = kw_ms1$p.value
)

############################################################
## Marker labels and order for panel C
############################################################
marker_label_map <- c(
  "CCL4"     = "CCL4",
  "CCL3"     = "CCL3",
  "CCL2"     = "CCL2",
  "IL_8_1"   = "CXCL8",
  "IL_1RA"   = "IL-1RA",
  "IL_10"    = "IL-10",
  "IL_6"     = "IL-6",
  "IL_1beta" = "IL-1\u03B2",
  "TNF"      = "TNF"
)

marker_order_display <- c(
  "CCL4",
  "CCL3",
  "CCL2",
  "IL_8_1",
  "IL_1RA",
  "IL_10",
  "IL_6",
  "IL_1beta",
  "TNF"
)

marker_levels_plot <- marker_order_display

############################################################
## Core function: Hedges' g + RAW P
############################################################
compute_hedges_lps <- function(
    studydata,
    wb_long,
    group_var,
    group_levels,
    axis_label
) {
  
  id_key <- studydata %>%
    transmute(
      EB_id = as.character(EB_id),
      axis_group = factor(.data[[group_var]], levels = group_levels)
    )
  
  wb <- wb_long %>%
    mutate(ID = as.character(ID)) %>%
    semi_join(id_key, by = c("ID" = "EB_id")) %>%
    left_join(id_key, by = c("ID" = "EB_id"))
  
  wb_clean <- wb %>%
    filter(
      stimulation == "LPS",
      !is.na(marker),
      !is.na(value),
      !is.na(axis_group)
    ) %>%
    mutate(
      value_raw = as.character(value),
      v_dot   = parse_number(value_raw, locale = locale(decimal_mark = ".", grouping_mark = ",")),
      v_comma = parse_number(value_raw, locale = locale(decimal_mark = ",", grouping_mark = ".")),
      value   = coalesce(v_dot, v_comma),
      value   = if_else(str_detect(value_raw, "^\\s*<"),
                        parse_number(str_replace(value_raw, "<","")),
                        value),
      value   = if_else(str_detect(value_raw, "^\\s*>"),
                        parse_number(str_replace(value_raw, ">","")),
                        value),
      value_num = log10(pmax(value, 1e-9)),
      axis_group = forcats::fct_relevel(axis_group, group_levels),
      marker = factor(marker, levels = marker_levels_plot)
    ) %>%
    dplyr::select(marker, value_num, axis_group)
  
  es <- wb_clean %>%
    group_by(marker) %>%
    rstatix::cohens_d(
      value_num ~ axis_group,
      var.equal = FALSE,
      hedges.correction = TRUE,
      ci = TRUE
    ) %>%
    as.data.frame() %>%
    ungroup()
  
  if ("effsize" %in% names(es)) {
    es <- dplyr::rename(es, g = effsize)
  } else if ("estimate" %in% names(es)) {
    es <- dplyr::rename(es, g = estimate)
  } else {
    stop("cohens_d() output has no 'effsize' or 'estimate'. Columns: ",
         paste(names(es), collapse = ", "))
  }
  
  tt <- wb_clean %>%
    group_by(marker) %>%
    rstatix::t_test(value_num ~ axis_group, var.equal = FALSE) %>%
    as.data.frame() %>%
    ungroup() %>%
    mutate(p_raw = p) %>%
    dplyr::select(marker, statistic, p_raw)
  
  es %>%
    left_join(tt, by = "marker") %>%
    mutate(
      axis = axis_label,
      p_BH = p_raw,
      p_BH.signif = case_when(
        is.na(p_raw) ~ "ns",
        p_raw < 1e-4 ~ "****",
        p_raw < 1e-3 ~ "***",
        p_raw < 1e-2 ~ "**",
        p_raw < 5e-2 ~ "*",
        TRUE ~ "ns"
      ),
      dir = case_when(
        statistic > 0 ~ "UP",
        statistic < 0 ~ "DOWN",
        TRUE ~ NA_character_
      ),
      signif_dir = ifelse(p_BH.signif == "ns" | is.na(dir), "ns", paste(p_BH.signif, dir)),
      marker = factor(as.character(marker), levels = marker_levels_plot)
    )
}

############################################################
## Compute Hedges' g within MS1_high for panel C
############################################################
studydata_ms1high <- studydata %>%
  filter(MS1_group == "MS1_high")

eff_ms1high_il1r2 <- compute_hedges_lps(
  studydata = studydata_ms1high,
  wb_long   = wb_long,
  group_var = "IL1R2_group",
  group_levels = c("IL1R2_positive", "IL1R2_zero"),
  axis_label = "IL1R2: positive vs zero (within MS1-high)"
) %>%
  mutate(
    axis   = factor(axis, levels = "IL1R2: positive vs zero (within MS1-high)"),
    marker = factor(marker, levels = marker_levels_plot)
  )

############################################################
## Forest helpers for panel C
############################################################
strip_bg_ms1high <- eff_ms1high_il1r2 %>%
  distinct(axis, marker) %>%
  arrange(axis, marker) %>%
  group_by(axis) %>%
  mutate(
    row_id = row_number(),
    ymin = row_id - 0.5,
    ymax = row_id + 0.5
  ) %>%
  filter(row_id %% 2 == 1) %>%
  ungroup()

color_mapping <- c(
  "ns"="#888888",
  "* UP"="#ffd4c4","** UP"="#FF9F94","*** UP"="#ff5c33","**** UP"="#cc2900",
  "* DOWN"="#c2d1ff","** DOWN"="#99b3ff","*** DOWN"="#668cff","**** DOWN"="#335CFF"
)

star_df_ms1high <- eff_ms1high_il1r2 %>%
  filter(p_BH.signif != "ns") %>%
  mutate(
    x_star = max(conf.high, na.rm = TRUE) +
      0.05 * diff(range(conf.low, conf.high, na.rm = TRUE))
  )

x_min_forest <- min(eff_ms1high_il1r2$conf.low,  na.rm = TRUE)
x_max_forest <- max(eff_ms1high_il1r2$conf.high, na.rm = TRUE)
x_pad_forest <- 0.24 * (x_max_forest - x_min_forest)

############################################################
## Panel C
############################################################
p_C <- ggplot(eff_ms1high_il1r2, aes(x = g, y = marker)) +
  geom_rect(
    data = strip_bg_ms1high,
    inherit.aes = FALSE,
    aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
    fill = "gray95"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.45) +
  geom_segment(aes(x = conf.low, xend = conf.high, yend = marker), linewidth = 0.45) +
  geom_point(aes(fill = signif_dir), shape = 22, size = 4.0, color = "black") +
  geom_text(
    data = star_df_ms1high,
    inherit.aes = FALSE,
    aes(x = x_star, y = marker, label = p_BH.signif),
    hjust = 0, size = STAR_SIZE,
    family = FONT_FAMILY, color = "black"
  ) +
  scale_fill_manual(values = color_mapping, na.value = "#888888") +
  scale_y_discrete(labels = marker_label_map) +
  labs(
    x = "Hedges’ g (95% CI)",
    y = NULL
  ) +
  coord_cartesian(
    xlim = c(x_min_forest - x_pad_forest, x_max_forest + x_pad_forest),
    clip = "off"
  ) +
  forest_theme +
  labs(tag = "C")

############################################################
## Reactome section for panel D
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

ms1_high_value <- "MS1_high"
il1r2_levels_within <- c("IL1R2_zero", "IL1R2_positive")

## IMPORTANT: Metabolism removed here
mother_keep <- c(
  "Innate Immune System",
  "Adaptive Immune System",
  "Cytokine Signaling in Immune system",
  "Hemostasis"
)

fdr_cutoff_union <- 0.05
bar_width <- 0.72

############################################################
## Helper: load counts as integer matrix
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
## Helper: align phenotype to counts
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
## Helper: top1000 signatures
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
## Helper: remove UNION signature genes
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
## Helper: run DESeq2 + Reactome GSEA
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
## Prepare Reactome backbone
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
## Load phenotype once
############################################################
studydata_rna <- readr::read_csv(pheno_file, show_col_types = FALSE) %>%
  dplyr::mutate(Patient_ID = as.character(Patient_ID))

############################################################
## Subset: MS1_high only, then align IL1R2 groups to counts
############################################################
study_ms1_high <- studydata_rna %>%
  dplyr::filter(MS1_group == ms1_high_value) %>%
  dplyr::filter(IL1R2_group %in% il1r2_levels_within) %>%
  dplyr::mutate(IL1R2_group = factor(IL1R2_group, levels = il1r2_levels_within))

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
## Circularity control
############################################################
cnt_within_cc <- remove_union_signature_top1000(
  cnt_int        = cnt_within,
  ms1_sig_file   = ms1_sig_file,
  il1r2_sig_file = il1r2_sig_file,
  n              = 1000
)

############################################################
## Run Reactome GSEA
############################################################
reactome_df_within <- run_deseq2_reactome(
  cnt_int      = cnt_within_cc,
  meta         = meta_within,
  contrast_vec = c("group", "IL1R2_positive", "IL1R2_zero")
)

############################################################
## Prepare GSEA table
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
## Child pathway selection
############################################################
all_children <- reactome_edges %>% dplyr::pull(Children) %>% unique()

master_children <- gsea_within %>%
  dplyr::filter(Description %in% all_children, !is.na(p.adjust), p.adjust < fdr_cutoff_union) %>%
  dplyr::pull(Description) %>%
  unique() %>%
  sort()

############################################################
## Build backbone rows
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
## Join GSEA results
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
## Panel D
############################################################
p_D <- ggplot(df_plot, aes(x = NES_plot, y = y_label)) +
  geom_vline(xintercept = 0, linewidth = 1.05, color = "black") +
  geom_col(
    aes(fill = bar_fill),
    width = bar_width,
    color = "black",
    linewidth = 0.25,
    na.rm = FALSE
  ) +
  scale_fill_gradient2(
    low  = "#4575b4",
    mid  = "white",
    high = "#d73027",
    limits = c(-nes_lim, nes_lim),
    name = "NES",
    na.value = "grey75"
  ) +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = NULL,
    tag = "D"
  ) +
  reactome_theme

############################################################
## Harmonize A and B + manual tags
############################################################
p_A <- p_A1 +
  labs(tag = "A") +
  theme(
    plot.margin = margin(5.5, 5.5, 4, 5.5, unit = "mm")
  )

p_B <- p_B1 +
  labs(tag = "B") +
  theme(
    plot.margin = margin(5.5, 5.5, 4, 5.5, unit = "mm")
  )

############################################################
## Build rows
############################################################
top_row <- (p_A | p_B) +
  plot_layout(widths = c(1, 1))

bottom_row <- (p_C | p_D) +
  plot_layout(widths = c(1.05, 1.55))

############################################################
## Final combined figure
############################################################
final_combined <- (top_row / free(bottom_row, side = "l")) +
  plot_layout(heights = c(1.00, 1.25)) &
  theme(
    text = element_text(family = FONT_FAMILY),
    plot.tag = element_text(
      family = FONT_FAMILY,
      size = TAG_SIZE,
      face = "plain"
    ),
    plot.tag.position = c(0.01, 0.98)
  )

print(final_combined)

############################################################
## Export
############################################################
ggsave(
  filename = "Figure/Figure_combined_final_clean_ABCD.svg",
  plot     = final_combined,
  width    = 330,
  height   = 205,
  units    = "mm",
  dpi      = 600,
  device   = svglite
)