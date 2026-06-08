## =========================================================
## Baseline plasma Hedges' g plots
## FINAL corrected version
## 1) Separate figures: markers ordered within each domain by MS1 or IL1R2 effect size
## 2) Combined figures: marker order is defined by MS1, then IL1R2 is matched to that order
## 3) No gap between domains
## 4) SVG text remains editable in Illustrator
## =========================================================

rm(list = ls())

## ===== Packages =====
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(rstatix)
  library(purrr)
  library(tibble)
  library(grid)
  library(svglite)
})

## =========================================================
## Global style (LOCKED)
## =========================================================
FONT_FAMILY <- "Arial"   # if Illustrator has trouble, try "Helvetica"
BASE_SIZE   <- 14
TITLE_SIZE  <- 14
TAG_SIZE    <- 16
AXIS_TITLE  <- 13
AXIS_TEXT   <- 11

LINE_AXIS   <- 0.6
LINE_TICKS  <- 0.6
MARGINS_MM  <- c(5.5, 5.5, 8.5, 5.5)

## =========================================================
## 0) Shared settings
## =========================================================

marker_label_map <- c(
  ## cytokines
  "IL-1a"         = "IL-1\u03b1",
  "IL-1b"         = "IL-1\u03b2",
  "TNFa"          = "TNF",
  "IFNa"          = "IFN\u03b1",
  "IFNg"          = "IFN\u03b3",
  
  ## chemokines
  "MIP-3a"        = "CCL20",
  "MIP-3b"        = "CCL19",
  "MIP-1a"        = "CCL3",
  "MIP-1b"        = "CCL4",
  "MCP-1"         = "CCL2",
  "IL-8"          = "CXCL8",
  "Rantes"        = "CCL5",
  "GROa"          = "CXCL1",
  "GROb"          = "CXCL2",
  "IP-10"         = "CXCL10",
  
  ## vascular / coagulation supplementary panel
  "VCAM"           = "VCAM-1",
  "E-Selectin"     = "E-selectin",
  "Angiopoetin-1"  = "Angiopoietin-1",
  "Angiopoietin-2" = "Angiopoietin-2",
  "D-Dimer"        = "D-dimer",
  "Tissue-Factor"  = "Tissue factor",
  "Thrombomodulin" = "Thrombomodulin"
)

make_display_label <- function(x) {
  ifelse(x %in% names(marker_label_map), marker_label_map[x], x)
}

## =========================================================
## Domains
## =========================================================

ordered_markers_main <- list(
  Cytokines = c(
    "IL-1a","IL-1b","IL-2","IL-3","IL-4","IL-5","IL-6","IL-7",
    "IL-13","IL-15","IL-17A","IL-18","IL-33","TNFa",
    "IFNa","IFNg","IL-10","IL-1RA","G-CSF","GM-CSF"
  ),
  Chemokines = c(
    "MCP-1","MIP-1a","MIP-1b","Rantes","MIP-3b","MIP-3a",
    "GROa","GROb","IL-8","IP-10"
  ),
  Inflammation = c(
    "Procalcitonin","CRP","Ferritin","TREM-1","CD163"
  )
)

ordered_markers_s8 <- list(
  "Vascular and coagulation responses" = c(
    "E-Selectin",
    "Angiopoetin-1",
    "Angiopoietin-2",
    "VCAM",
    "Syndecan-1",
    "Thrombomodulin",
    "D-Dimer",
    "Tissue-Factor"
  )
)

domain_map_main <- ordered_markers_main
domain_map_s8   <- ordered_markers_s8

domain_order_main <- c("Cytokines", "Chemokines", "Inflammation")
domain_order_s8   <- c("Vascular and coagulation responses")

color_mapping <- c(
  "ns"        = "#B0B0B0",
  "* UP"      = "#FFD4C4",
  "** UP"     = "#FF9F94",
  "*** UP"    = "#FF5C33",
  "**** UP"   = "#CC2900",
  "* DOWN"    = "#C2D1FF",
  "** DOWN"   = "#99B3FF",
  "*** DOWN"  = "#668CFF",
  "**** DOWN" = "#335CFF"
)

## =========================================================
## 1) Helper functions
## =========================================================

prepare_bio_long <- function(bio_df, domain_map) {
  domain_df <- enframe(domain_map, name = "domain", value = "marker") %>%
    unnest(marker)
  
  bio_df %>%
    pivot_longer(
      cols = where(is.numeric),
      names_to = "marker",
      values_to = "value"
    ) %>%
    filter(!is.na(value)) %>%
    mutate(value_num = log10(pmax(value, 1e-9))) %>%
    inner_join(domain_df, by = "marker")
}

calc_effect_df <- function(bio_long, group_var, domain_levels) {
  bio_long %>%
    group_by(domain, marker) %>%
    rstatix::cohens_d(
      formula = as.formula(paste0("value_num ~ ", group_var)),
      var.equal = FALSE,
      hedges.correction = TRUE,
      ci = TRUE
    ) %>%
    ungroup() %>%
    dplyr::rename(g = effsize) %>%
    left_join(
      bio_long %>%
        group_by(domain, marker) %>%
        rstatix::t_test(
          formula = as.formula(paste0("value_num ~ ", group_var)),
          var.equal = FALSE
        ) %>%
        ungroup() %>%
        mutate(p_BH = p.adjust(p, method = "BH")) %>%
        dplyr::select(domain, marker, p_BH),
      by = c("domain", "marker")
    ) %>%
    mutate(
      signif = case_when(
        p_BH < 1e-4 ~ "****",
        p_BH < 1e-3 ~ "***",
        p_BH < 1e-2 ~ "**",
        p_BH < 0.05 ~ "*",
        TRUE ~ ""
      ),
      domain = factor(domain, levels = domain_levels),
      marker_display = make_display_label(marker)
    )
}

## ---------------------------------------------------------
## Order for separate plots: sort by each figure's own effect size
## ---------------------------------------------------------
add_effect_order_separate <- function(eff_df) {
  out <- eff_df %>%
    group_by(domain) %>%
    arrange(desc(g), marker, .by_group = TRUE) %>%
    mutate(
      row_id = row_number(),
      n_rows = n(),
      plot_pos = n_rows - row_id + 1,
      bg_flag = row_id %% 2,
      marker_id = paste0(as.character(domain), "___", marker)
    ) %>%
    ungroup()
  
  marker_levels_topdown <- out %>%
    arrange(domain, row_id) %>%
    pull(marker_id) %>%
    unique()
  
  out %>%
    mutate(marker_id = factor(marker_id, levels = rev(marker_levels_topdown)))
}

## ---------------------------------------------------------
## Create master order from MS1 for combined plots
## ---------------------------------------------------------
make_master_order_from_ms1 <- function(eff_df_ms1) {
  eff_df_ms1 %>%
    group_by(domain) %>%
    arrange(desc(g), marker, .by_group = TRUE) %>%
    mutate(row_id = row_number()) %>%
    ungroup() %>%
    dplyr::select(domain, marker, row_id)
}

## ---------------------------------------------------------
## Apply external master order (from MS1) to any effect table
## ---------------------------------------------------------
apply_master_order <- function(eff_df, order_df) {
  out <- eff_df %>%
    left_join(order_df, by = c("domain", "marker")) %>%
    group_by(domain) %>%
    arrange(row_id, .by_group = TRUE) %>%
    mutate(
      n_rows = n(),
      plot_pos = n_rows - row_id + 1,
      bg_flag = row_id %% 2,
      marker_id = paste0(as.character(domain), "___", marker)
    ) %>%
    ungroup()
  
  marker_levels_topdown <- order_df %>%
    arrange(domain, row_id) %>%
    mutate(marker_id = paste0(as.character(domain), "___", marker)) %>%
    pull(marker_id) %>%
    unique()
  
  out %>%
    mutate(marker_id = factor(marker_id, levels = rev(marker_levels_topdown)))
}

base_forest_theme <- function(show_legend = FALSE) {
  theme_classic(base_size = BASE_SIZE) +
    theme(
      text = element_text(family = FONT_FAMILY),
      axis.title.x = element_text(size = AXIS_TITLE),
      axis.text.x  = element_text(size = AXIS_TEXT),
      axis.text.y  = element_text(size = AXIS_TEXT, face = "bold"),
      strip.text.x = element_text(size = TAG_SIZE, face = "plain"),
      strip.text.y = element_text(size = TAG_SIZE, face = "plain", angle = 270),
      strip.background = element_blank(),
      axis.line.x = element_line(linewidth = LINE_AXIS),
      axis.line.y = element_blank(),
      axis.ticks  = element_line(linewidth = LINE_TICKS),
      axis.ticks.y = element_blank(),
      panel.spacing.x = unit(1.0, "lines"),
      panel.spacing.y = unit(0, "lines"),
      legend.position = if (show_legend) "top" else "none",
      legend.text = element_text(size = AXIS_TEXT),
      plot.margin = margin(
        t = MARGINS_MM[1], r = MARGINS_MM[2],
        b = MARGINS_MM[3], l = MARGINS_MM[4],
        unit = "mm"
      )
    )
}

make_single_plot <- function(eff_df,
                             file_out,
                             width = 8,
                             height = 8,
                             show_legend = FALSE) {
  
  strip_bg <- eff_df %>%
    filter(bg_flag == 1) %>%
    distinct(domain, marker_id, plot_pos) %>%
    mutate(
      ymin = plot_pos - 0.5,
      ymax = plot_pos + 0.5
    )
  
  star_df <- eff_df %>%
    filter(signif != "") %>%
    group_by(domain) %>%
    mutate(
      x_lab = max(conf.high, na.rm = TRUE) + 0.08
    ) %>%
    ungroup()
  
  p <- ggplot(eff_df, aes(x = g, y = marker_id)) +
    geom_rect(
      data = strip_bg,
      inherit.aes = FALSE,
      aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
      fill = "grey90"
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
    geom_segment(
      aes(x = conf.low, xend = conf.high, yend = marker_id),
      linewidth = 0.4
    ) +
    geom_point(
      aes(fill = sig_dir),
      shape = 22,
      size = 4.2,
      color = "black"
    ) +
    geom_text(
      data = star_df,
      inherit.aes = FALSE,
      aes(x = x_lab, y = marker_id, label = signif),
      hjust = 0,
      size = AXIS_TEXT / ggplot2::.pt,
      family = FONT_FAMILY,
      color = "black"
    ) +
    scale_fill_manual(
      values = color_mapping,
      drop = FALSE,
      guide = if (show_legend) guide_legend(title = NULL) else "none"
    ) +
    scale_y_discrete(
      labels = function(x) {
        raw_marker <- sub("^.*___", "", x)
        make_display_label(raw_marker)
      },
      expand = expansion(add = c(0.5, 0.5))
    ) +
    facet_grid(domain ~ ., scales = "free_y", space = "free_y") +
    labs(
      x = "Hedges’ g, 95% CI",
      y = NULL
    ) +
    coord_cartesian(clip = "off") +
    base_forest_theme(show_legend = show_legend)
  
  print(p)
  
  ggsave(
    filename = file_out,
    plot     = p,
    device   = svglite::svglite,
    width    = width,
    height   = height
  )
  
  invisible(p)
}

make_combined_plot <- function(plot_df,
                               file_out,
                               width = 10,
                               height = 8) {
  
  strip_bg <- plot_df %>%
    filter(bg_flag == 1) %>%
    distinct(domain, marker_id, plot_pos) %>%
    mutate(
      ymin = plot_pos - 0.5,
      ymax = plot_pos + 0.5
    )
  
  star_df <- plot_df %>%
    filter(signif != "") %>%
    group_by(domain, contrast) %>%
    mutate(
      x_lab = max(conf.high, na.rm = TRUE) + 0.08
    ) %>%
    ungroup()
  
  p <- ggplot(plot_df, aes(x = g, y = marker_id)) +
    geom_rect(
      data = strip_bg,
      inherit.aes = FALSE,
      aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
      fill = "grey90"
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
    geom_segment(
      aes(x = conf.low, xend = conf.high, yend = marker_id),
      linewidth = 0.4
    ) +
    geom_point(
      aes(fill = sig_dir),
      shape = 22,
      size = 4.2,
      color = "black"
    ) +
    geom_text(
      data = star_df,
      inherit.aes = FALSE,
      aes(x = x_lab, y = marker_id, label = signif),
      hjust = 0,
      size = AXIS_TEXT / ggplot2::.pt,
      family = FONT_FAMILY,
      color = "black"
    ) +
    facet_grid(
      domain ~ contrast,
      scales = "free_y",
      space  = "free_y",
      labeller = labeller(
        contrast = c(
          "MS1_high_vs_low"        = "A   MS1-high",
          "IL1R2_positive_vs_zero" = "B   IL1R2-positive"
        )
      )
    ) +
    scale_fill_manual(
      values = color_mapping,
      drop = FALSE,
      guide = "none"
    ) +
    scale_y_discrete(
      labels = function(x) {
        raw_marker <- sub("^.*___", "", x)
        make_display_label(raw_marker)
      },
      expand = expansion(add = c(0.5, 0.5))
    ) +
    labs(
      x = "Hedges’ g, 95% CI",
      y = NULL
    ) +
    coord_cartesian(clip = "off") +
    base_forest_theme(show_legend = FALSE)
  
  print(p)
  
  ggsave(
    filename = file_out,
    plot     = p,
    device   = svglite::svglite,
    width    = width,
    height   = height
  )
  
  invisible(p)
}

## =========================================================
## 2) Load data
## =========================================================

studydata_with_ms1 <- read_csv(
  "Original_data/studydata_with_ms1.csv",
  show_col_types = FALSE
)

studydata_with_il1r2 <- read_csv(
  "Original_data/studydata_with_il1r2.csv",
  show_col_types = FALSE
)

bio_raw <- read.table(
  "Original_Data/EB_optimact_final_biomarkers_luminex_cba.csv",
  header = TRUE,
  sep = ";",
  dec = ",",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

bio_raw <- bio_raw %>%
  filter(str_starts(id, "EB_")) %>%
  mutate(ID = str_remove(id, "^EB_"))

## =========================================================
## 3) Build datasets
## =========================================================

bio_ms1 <- bio_raw %>%
  left_join(
    studydata_with_ms1 %>%
      transmute(
        ID = as.character(EB_id),
        MS1_group = factor(MS1_group, levels = c("MS1_high", "MS1_low"))
      ),
    by = "ID"
  ) %>%
  filter(!is.na(MS1_group))

bio_il1r2 <- bio_raw %>%
  left_join(
    studydata_with_il1r2 %>%
      transmute(
        ID = as.character(EB_id),
        IL1R2_group = factor(
          IL1R2_group,
          levels = c("IL1R2_positive", "IL1R2_zero")
        )
      ),
    by = "ID"
  ) %>%
  filter(!is.na(IL1R2_group))

bio_long_ms1_main   <- prepare_bio_long(bio_ms1,   domain_map_main)
bio_long_il1r2_main <- prepare_bio_long(bio_il1r2, domain_map_main)

bio_long_ms1_s8     <- prepare_bio_long(bio_ms1,   domain_map_s8)
bio_long_il1r2_s8   <- prepare_bio_long(bio_il1r2, domain_map_s8)

## =========================================================
## 4) Calculate effect tables
## =========================================================

eff_ms1_main_raw <- calc_effect_df(
  bio_long      = bio_long_ms1_main,
  group_var     = "MS1_group",
  domain_levels = domain_order_main
) %>%
  mutate(
    contrast = "MS1_high_vs_low",
    sig_dir = case_when(
      signif == "" ~ "ns",
      g > 0 ~ paste0(signif, " UP"),
      g < 0 ~ paste0(signif, " DOWN")
    )
  )

eff_il1r2_main_raw <- calc_effect_df(
  bio_long      = bio_long_il1r2_main,
  group_var     = "IL1R2_group",
  domain_levels = domain_order_main
) %>%
  mutate(
    contrast = "IL1R2_positive_vs_zero",
    sig_dir = case_when(
      signif == "" ~ "ns",
      g > 0 ~ paste0(signif, " UP"),
      g < 0 ~ paste0(signif, " DOWN")
    )
  )

eff_ms1_s8_raw <- calc_effect_df(
  bio_long      = bio_long_ms1_s8,
  group_var     = "MS1_group",
  domain_levels = domain_order_s8
) %>%
  mutate(
    contrast = "MS1_high_vs_low",
    sig_dir = case_when(
      signif == "" ~ "ns",
      g > 0 ~ paste0(signif, " UP"),
      g < 0 ~ paste0(signif, " DOWN")
    )
  )

eff_il1r2_s8_raw <- calc_effect_df(
  bio_long      = bio_long_il1r2_s8,
  group_var     = "IL1R2_group",
  domain_levels = domain_order_s8
) %>%
  mutate(
    contrast = "IL1R2_positive_vs_zero",
    sig_dir = case_when(
      signif == "" ~ "ns",
      g > 0 ~ paste0(signif, " UP"),
      g < 0 ~ paste0(signif, " DOWN")
    )
  )

## =========================================================
## 5) Separate figures
## each ordered by its own effect size
## =========================================================

eff_ms1_main <- add_effect_order_separate(eff_ms1_main_raw)
eff_il1r2_main <- add_effect_order_separate(eff_il1r2_main_raw)

eff_ms1_s8 <- add_effect_order_separate(eff_ms1_s8_raw)
eff_il1r2_s8 <- add_effect_order_separate(eff_il1r2_s8_raw)

p_ms1_main <- make_single_plot(
  eff_df   = eff_ms1_main,
  file_out = "Figure/baseline_plasma_MAIN_MS1_high_vs_low.svg",
  width    = 8,
  height   = 8,
  show_legend = FALSE
)

p_il1r2_main <- make_single_plot(
  eff_df   = eff_il1r2_main,
  file_out = "Figure/baseline_plasma_MAIN_IL1R2_positive_vs_zero.svg",
  width    = 8,
  height   = 8,
  show_legend = FALSE
)

## =========================================================
## 6) Combined MAIN figure
## order defined by MS1, IL1R2 matched to same order
## =========================================================

master_order_main <- make_master_order_from_ms1(eff_ms1_main_raw)

eff_ms1_main_combined <- apply_master_order(eff_ms1_main_raw, master_order_main)
eff_il1r2_main_combined <- apply_master_order(eff_il1r2_main_raw, master_order_main)

plot_df_main <- bind_rows(
  eff_ms1_main_combined,
  eff_il1r2_main_combined
) %>%
  mutate(
    contrast = factor(
      contrast,
      levels = c("MS1_high_vs_low", "IL1R2_positive_vs_zero")
    ),
    domain = factor(domain, levels = domain_order_main)
  )

p_main_combined <- make_combined_plot(
  plot_df  = plot_df_main,
  file_out = "Figure/Figure4_baseline_plasma_MAIN_MS1_IL1R2_combined.svg",
  width    = 10,
  height   = 8.5
)

## =========================================================
## 7) Combined S8 figure
## order defined by MS1, IL1R2 matched to same order
## =========================================================

master_order_s8 <- make_master_order_from_ms1(eff_ms1_s8_raw)

eff_ms1_s8_combined <- apply_master_order(eff_ms1_s8_raw, master_order_s8)
eff_il1r2_s8_combined <- apply_master_order(eff_il1r2_s8_raw, master_order_s8)

plot_df_s8 <- bind_rows(
  eff_ms1_s8_combined,
  eff_il1r2_s8_combined
) %>%
  mutate(
    contrast = factor(
      contrast,
      levels = c("MS1_high_vs_low", "IL1R2_positive_vs_zero")
    ),
    domain = factor(domain, levels = domain_order_s8)
  )

p_s8_combined <- make_combined_plot(
  plot_df  = plot_df_s8,
  file_out = "Figure/FigureS8_vascular_coagulation_MS1_IL1R2_combined.svg",
  width    = 10,
  height   = 4.5
)

## =========================================================
## 8) Save tables if needed
## =========================================================

write_csv(
  eff_ms1_main_raw %>%
    mutate(marker_display = make_display_label(marker)) %>%
    dplyr::select(domain, marker, marker_display, g, conf.low, conf.high, p_BH, signif),
  "Figure/table_MAIN_MS1_effect_sizes.csv"
)

write_csv(
  eff_il1r2_main_raw %>%
    mutate(marker_display = make_display_label(marker)) %>%
    dplyr::select(domain, marker, marker_display, g, conf.low, conf.high, p_BH, signif),
  "Figure/table_MAIN_IL1R2_effect_sizes.csv"
)

write_csv(
  eff_ms1_s8_raw %>%
    mutate(marker_display = make_display_label(marker)) %>%
    dplyr::select(domain, marker, marker_display, g, conf.low, conf.high, p_BH, signif),
  "Figure/table_S8_MS1_effect_sizes.csv"
)

write_csv(
  eff_il1r2_s8_raw %>%
    mutate(marker_display = make_display_label(marker)) %>%
    dplyr::select(domain, marker, marker_display, g, conf.low, conf.high, p_BH, signif),
  "Figure/table_S8_IL1R2_effect_sizes.csv"
)

## =========================================================
## END
## =========================================================