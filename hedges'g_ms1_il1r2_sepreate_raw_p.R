###########################################################
## Figure:
## a/b) Forest plot as one faceted panel: MS1 and IL1R2+
## c) LPS-induced TNF vs percentage of MS1 cells
## d) LPS-induced TNF across IL1R2 strata
##
## ERJ-oriented settings:
## - figure size: 180 mm x 220 mm
## - width kept at ERJ maximum: 180 mm
## - top forest plot spans full width
## - bottom row has c) and d)
## - panel labels NOT bold
## - global text size matched to distribution figure
## - much wider spacing between the two forest facets
############################################################

rm(list = ls())

############################################################
## Packages
############################################################
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(rstatix)
  library(mgcv)
  library(patchwork)
  library(tidyr)
  library(grid)
  library(svglite)
})

############################################################
## Input
############################################################
studydata <- read_csv(
  "Original_data/studydata_with_ms1_il1r2.csv",
  show_col_types = FALSE
)

wb_long <- read_csv(
  "Original_data/whole_blood_stimuli_long.csv",
  show_col_types = FALSE
)

if (!dir.exists("Figure")) dir.create("Figure")

############################################################
## ERJ-safe export size
############################################################
FIG_WIDTH_MM  <- 180
FIG_HEIGHT_MM <- 220

stopifnot(FIG_WIDTH_MM  <= 180)
stopifnot(FIG_HEIGHT_MM <= 230)

############################################################
## Global style
## Same as final distribution figure
############################################################
FONT_FAMILY <- "Arial"
BASE_SIZE   <- 14
TITLE_SIZE  <- 14
TAG_SIZE    <- 16
AXIS_TITLE  <- 13
AXIS_TEXT   <- 11

LINE_AXIS   <- 0.6
LINE_TICKS  <- 0.6
POINT_SIZE  <- 2.0
BIN_EDGE_W  <- 0.2
MARGINS_MM  <- c(5.5, 5.5, 5.5, 5.5)

## Slightly smaller only for top-arrow text to avoid crowding
ARROW_TEXT_SIZE  <- 8.5
DOMAIN_TEXT_SIZE <- 8.5
STAR_SIZE        <- 3.2

## Forest spacing
## Previous value was 1.6 lines; this is twice as wide.
FOREST_FACET_SPACING <- 3.2

col_ms1   <- "#2F4B7C"
col_il1r2 <- "#A23B3B"

base_theme <- theme_classic(base_size = BASE_SIZE) +
  theme(
    text        = element_text(family = FONT_FAMILY),
    plot.title  = element_text(size = TITLE_SIZE, face = "plain", hjust = 0),
    plot.tag    = element_text(size = TAG_SIZE, face = "plain", family = FONT_FAMILY),
    axis.title  = element_text(size = AXIS_TITLE),
    axis.text   = element_text(size = AXIS_TEXT),
    axis.line   = element_line(linewidth = LINE_AXIS),
    axis.ticks  = element_line(linewidth = LINE_TICKS),
    plot.margin = margin(
      MARGINS_MM[1],
      MARGINS_MM[2],
      MARGINS_MM[3],
      MARGINS_MM[4],
      unit = "mm"
    )
  )

############################################################
## Marker labels and order
############################################################
marker_order_raw <- c(
  "TNF",
  "IL_1beta",
  "IL_6",
  "IL_10",
  "IL_1RA",
  "IL_8_1",
  "CCL2",
  "CCL3",
  "CCL4"
)

marker_label_map <- c(
  "TNF"      = "TNF",
  "IL_1beta" = "IL-1\u03B2",
  "IL_6"     = "IL-6",
  "IL_10"    = "IL-10",
  "IL_1RA"   = "IL-1RA",
  "IL_8_1"   = "CXCL8",
  "CCL2"     = "CCL2",
  "CCL3"     = "CCL3",
  "CCL4"     = "CCL4"
)

############################################################
## Helper: robust parsing
############################################################
parse_wb_value <- function(x) {
  x_chr <- as.character(x)
  
  v_dot <- parse_number(
    x_chr,
    locale = locale(decimal_mark = ".", grouping_mark = ",")
  )
  
  v_comma <- parse_number(
    x_chr,
    locale = locale(decimal_mark = ",", grouping_mark = ".")
  )
  
  value_num <- coalesce(v_dot, v_comma)
  
  value_num <- if_else(
    str_detect(x_chr, "^\\s*<"),
    parse_number(str_replace(x_chr, "<", "")),
    value_num
  )
  
  value_num <- if_else(
    str_detect(x_chr, "^\\s*>"),
    parse_number(str_replace(x_chr, ">", "")),
    value_num
  )
  
  value_num
}

############################################################
## Helper: Hedges' g + raw p
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
      marker %in% marker_order_raw,
      !is.na(marker),
      !is.na(value),
      !is.na(axis_group)
    ) %>%
    mutate(
      value_num   = parse_wb_value(value),
      value_log10 = log10(pmax(value_num, 1e-9)),
      axis_group  = forcats::fct_relevel(axis_group, group_levels),
      marker      = factor(marker, levels = marker_order_raw)
    ) %>%
    select(marker, value_log10, axis_group)
  
  es <- wb_clean %>%
    group_by(marker) %>%
    rstatix::cohens_d(
      value_log10 ~ axis_group,
      var.equal = FALSE,
      hedges.correction = TRUE,
      ci = TRUE
    ) %>%
    as.data.frame() %>%
    ungroup()
  
  if ("effsize" %in% names(es)) {
    es <- rename(es, g = effsize)
  } else if ("estimate" %in% names(es)) {
    es <- rename(es, g = estimate)
  } else {
    stop("cohens_d() output has no 'effsize' or 'estimate'")
  }
  
  tt <- wb_clean %>%
    group_by(marker) %>%
    rstatix::t_test(value_log10 ~ axis_group, var.equal = FALSE) %>%
    as.data.frame() %>%
    ungroup() %>%
    mutate(p_raw = p) %>%
    select(marker, statistic, p_raw)
  
  es %>%
    left_join(tt, by = "marker") %>%
    mutate(
      axis = axis_label,
      p_signif = case_when(
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
      signif_dir = ifelse(
        p_signif == "ns" | is.na(dir),
        "ns",
        paste(p_signif, dir)
      ),
      marker = factor(as.character(marker), levels = marker_order_raw)
    )
}

############################################################
## Forest data
############################################################
eff_all <- bind_rows(
  compute_hedges_lps(
    studydata = studydata,
    wb_long   = wb_long,
    group_var = "MS1_group",
    group_levels = c("MS1_high", "MS1_low"),
    axis_label   = "MS1"
  ),
  compute_hedges_lps(
    studydata = studydata,
    wb_long   = wb_long,
    group_var = "IL1R2_group",
    group_levels = c("IL1R2_positive", "IL1R2_zero"),
    axis_label   = "IL1R2"
  )
) %>%
  mutate(
    axis = factor(axis, levels = c("MS1", "IL1R2")),
    marker = factor(marker, levels = marker_order_raw)
  )

############################################################
## Forest helper data
############################################################
strip_bg <- eff_all %>%
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
  "* UP"="#ffd4c4",
  "** UP"="#FF9F94",
  "*** UP"="#ff5c33",
  "**** UP"="#cc2900",
  "* DOWN"="#c2d1ff",
  "** DOWN"="#99b3ff",
  "*** DOWN"="#668cff",
  "**** DOWN"="#335CFF"
)

x_rng <- eff_all %>%
  group_by(axis) %>%
  summarise(
    x_min = min(conf.low, na.rm = TRUE),
    x_max = max(conf.high, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    x_pad = 0.20 * (x_max - x_min),
    x_left = x_min - x_pad,
    x_right = x_max + x_pad
  )

## Stars kept inside each facet
star_df <- eff_all %>%
  filter(p_signif != "ns") %>%
  left_join(x_rng, by = "axis") %>%
  mutate(
    x_star = x_right - 0.04 * (x_right - x_left)
  )

arrow_bar_df <- x_rng %>%
  transmute(
    axis,
    xstart = x_left + 0.02 * (x_right - x_left),
    xend   = x_right - 0.02 * (x_right - x_left),
    y      = Inf
  )

arrow_text_df <- x_rng %>%
  mutate(
    x_left_mid  = (x_left + 0) / 2,
    x_right_mid = (0 + x_right) / 2
  ) %>%
  transmute(
    axis,
    x_left_mid,
    x_right_mid,
    y = Inf,
    lab_left = case_when(
      axis == "MS1"   ~ "Lower in MS1-high",
      axis == "IL1R2" ~ "Lower in IL1R2+",
      TRUE ~ "Lower"
    ),
    lab_right = case_when(
      axis == "MS1"   ~ "Higher in MS1-high",
      axis == "IL1R2" ~ "Higher in IL1R2+",
      TRUE ~ "Higher"
    )
  )

domain_text_df <- x_rng %>%
  transmute(
    axis,
    x = x_left + 0.03 * (x_right - x_left)
  ) %>%
  tidyr::crossing(
    tibble(
      marker = c("TNF", "IL_8_1"),
      domain_label = c("Cytokines", "Chemokines")
    )
  ) %>%
  mutate(
    marker = factor(marker, levels = marker_order_raw)
  )

############################################################
## Top forest plot: one faceted plot spanning full width
## This keeps marker names only once on the left
############################################################
p_forest <- ggplot(eff_all, aes(x = g, y = marker)) +
  geom_rect(
    data = strip_bg,
    inherit.aes = FALSE,
    aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
    fill = "gray95"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35) +
  geom_segment(
    data = arrow_bar_df,
    inherit.aes = FALSE,
    aes(x = xstart, xend = xend, y = y, yend = y),
    linewidth = 0.45,
    lineend = "butt",
    arrow = arrow(
      ends   = "both",
      type   = "closed",
      angle  = 25,
      length = unit(1.8, "mm")
    )
  ) +
  geom_text(
    data = arrow_text_df,
    inherit.aes = FALSE,
    aes(x = x_left_mid, y = y, label = lab_left),
    hjust = 0.5,
    vjust = 1.55,
    size = ARROW_TEXT_SIZE / ggplot2::.pt,
    family = FONT_FAMILY,
    fontface = "plain"
  ) +
  geom_text(
    data = arrow_text_df,
    inherit.aes = FALSE,
    aes(x = x_right_mid, y = y, label = lab_right),
    hjust = 0.5,
    vjust = 1.55,
    size = ARROW_TEXT_SIZE / ggplot2::.pt,
    family = FONT_FAMILY,
    fontface = "plain"
  ) +
  geom_text(
    data = domain_text_df,
    inherit.aes = FALSE,
    aes(x = x, y = marker, label = domain_label),
    hjust = 0,
    vjust = 1.65,
    size = DOMAIN_TEXT_SIZE / ggplot2::.pt,
    family = FONT_FAMILY,
    fontface = "plain"
  ) +
  geom_segment(
    aes(x = conf.low, xend = conf.high, yend = marker),
    linewidth = 0.35
  ) +
  geom_point(
    aes(fill = signif_dir),
    shape = 22,
    size = 3.4,
    color = "black",
    linewidth = 0.25
  ) +
  geom_text(
    data = star_df,
    inherit.aes = FALSE,
    aes(x = x_star, y = marker, label = p_signif),
    hjust = 0,
    size = STAR_SIZE,
    family = FONT_FAMILY,
    color = "black"
  ) +
  scale_fill_manual(values = color_mapping, na.value = "#888888") +
  scale_y_discrete(
    limits = rev(marker_order_raw),
    labels = marker_label_map
  ) +
  facet_grid(cols = vars(axis), scales = "free_x") +
  labs(
    tag = "a)",
    x = "Hedges’ g (95% CI)",
    y = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = BASE_SIZE) +
  theme(
    text = element_text(family = FONT_FAMILY),
    
    plot.tag = element_text(
      family = FONT_FAMILY,
      size = TAG_SIZE,
      face = "plain"
    ),
    plot.tag.position = c(0.005, 1.02),
    
    axis.text.y = element_text(size = AXIS_TEXT, face = "plain"),
    axis.text.x = element_text(size = AXIS_TEXT),
    axis.title.x = element_text(size = AXIS_TITLE),
    
    axis.line.x = element_line(linewidth = LINE_AXIS),
    axis.ticks  = element_line(linewidth = LINE_TICKS),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.position = "none",
    
    ## Much wider distance between the two forest facets
    panel.spacing.x = unit(FOREST_FACET_SPACING, "lines"),
    
    plot.margin = margin(
      t = 15,
      r = 3,
      b = 6,
      l = 3,
      unit = "mm"
    )
  )

############################################################
## Prepare TNF data
############################################################
df_tnf <- wb_long %>%
  filter(stimulation == "LPS", marker == "TNF") %>%
  mutate(
    ID = as.character(ID),
    value_num = parse_wb_value(value)
  ) %>%
  left_join(
    studydata %>%
      transmute(
        EB_id = as.character(EB_id),
        MS1_percent = as.numeric(MS1_percent),
        IL1R2_percent = as.numeric(IL1R2_percent)
      ),
    by = c("ID" = "EB_id")
  ) %>%
  filter(
    !is.na(value_num),
    value_num > 0,
    !is.na(MS1_percent),
    !is.na(IL1R2_percent)
  ) %>%
  mutate(
    value_log10 = log10(value_num)
  )

############################################################
## Panel c: TNF vs MS1
############################################################
p_c <- ggplot(df_tnf, aes(x = MS1_percent, y = value_log10)) +
  geom_point(color = col_ms1, alpha = 0.65, size = POINT_SIZE) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, k = 4),
    color = col_ms1,
    fill = col_ms1,
    alpha = 0.22,
    linewidth = 0.9
  ) +
  labs(
    tag = "c)",
    x = "Percentage of MS1 cells",
    y = "LPS-induced TNF (log10-transformed pg/mL)"
  ) +
  base_theme +
  theme(
    plot.tag.position = c(0.005, 1.02)
  )

############################################################
## Build IL1R2 strata for panel d
## 0 = IL1R2-zero
## 1-6 = strata among IL1R2+ patients
############################################################
df_il1r2_grp <- df_tnf %>%
  mutate(IL1R2_pos = IL1R2_percent > 0)

pos_breaks <- quantile(
  df_il1r2_grp$IL1R2_percent[df_il1r2_grp$IL1R2_percent > 0],
  probs = seq(0, 1, length.out = 7),
  na.rm = TRUE
)

pos_breaks <- unique(as.numeric(pos_breaks))

if (length(pos_breaks) < 7) {
  pos_breaks <- pretty(
    range(
      df_il1r2_grp$IL1R2_percent[df_il1r2_grp$IL1R2_percent > 0],
      na.rm = TRUE
    ),
    n = 6
  )
}

df_il1r2_grp <- df_il1r2_grp %>%
  mutate(
    IL1R2_group7 = case_when(
      IL1R2_percent == 0 ~ 0L,
      IL1R2_percent > 0  ~ as.integer(cut(
        IL1R2_percent,
        breaks = pos_breaks,
        include.lowest = TRUE,
        labels = 1:6
      )),
      TRUE ~ NA_integer_
    ),
    IL1R2_group7 = factor(IL1R2_group7, levels = 0:6),
    IL1R2_group7_num = as.numeric(as.character(IL1R2_group7))
  ) %>%
  filter(!is.na(IL1R2_group7_num))

############################################################
## GAM curve for panel d
## Fitted only over positive strata 1-6
############################################################
df_for_curve <- df_il1r2_grp %>%
  filter(IL1R2_group7_num >= 1)

gam_grp <- mgcv::gam(
  value_log10 ~ s(IL1R2_group7_num, k = 4),
  data = df_for_curve,
  method = "REML"
)

new_curve <- data.frame(
  IL1R2_group7_num = seq(1, 6, length.out = 200)
)

pred <- predict(gam_grp, newdata = new_curve, se.fit = TRUE)

curve_df <- new_curve %>%
  mutate(
    fit  = pred$fit,
    low  = pred$fit - 1.96 * pred$se.fit,
    high = pred$fit + 1.96 * pred$se.fit
  )

############################################################
## Panel d: TNF across IL1R2 strata
############################################################
p_d <- ggplot(df_il1r2_grp, aes(x = IL1R2_group7_num, y = value_log10)) +
  geom_boxplot(
    aes(group = IL1R2_group7_num),
    outlier.shape = NA,
    linewidth = 0.65,
    fill = "#F3D7CF"
  ) +
  geom_jitter(
    width = 0.18,
    height = 0,
    alpha = 0.55,
    size = POINT_SIZE,
    color = col_il1r2
  ) +
  geom_ribbon(
    data = curve_df,
    inherit.aes = FALSE,
    aes(x = IL1R2_group7_num, ymin = low, ymax = high),
    fill = col_il1r2,
    alpha = 0.20
  ) +
  geom_line(
    data = curve_df,
    inherit.aes = FALSE,
    aes(x = IL1R2_group7_num, y = fit),
    color = col_il1r2,
    linewidth = 0.9
  ) +
  scale_x_continuous(
    breaks = 0:6,
    labels = 0:6,
    limits = c(-0.2, 6.2)
  ) +
  labs(
    tag = "d)",
    x = "IL1R2 strata",
    y = NULL
  ) +
  base_theme +
  theme(
    plot.tag.position = c(0.005, 1.02)
  )

############################################################
## Combine figure
############################################################
bottom_row <- p_c + p_d +
  plot_layout(
    ncol = 2,
    widths = c(1, 1)
  )

final_plot <- p_forest / bottom_row +
  plot_layout(
    heights = c(1.55, 1.00)
  )

print(final_plot)

############################################################
## Export
############################################################
ggsave(
  filename = "Figure/Figure_LPS_MS1_IL1R2_combined_4panel.svg",
  plot     = final_plot,
  device   = svglite,
  width    = FIG_WIDTH_MM,
  height   = FIG_HEIGHT_MM,
  units    = "mm"
)

ggsave(
  filename = "Figure/Figure_LPS_MS1_IL1R2_combined_4panel.pdf",
  plot     = final_plot,
  width    = FIG_WIDTH_MM,
  height   = FIG_HEIGHT_MM,
  units    = "mm",
  device   = cairo_pdf
)