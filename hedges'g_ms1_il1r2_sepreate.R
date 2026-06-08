############################################################
## 0) Clean workspace
############################################################
rm(list = ls())

############################################################
## 1) Packages
############################################################
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(rstatix)
  library(purrr)
  library(grid)      # unit(), arrow()
  library(mgcv)      # gam
  library(patchwork) # combine panels
})

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
MARGINS_MM  <- c(5.5, 5.5, 8.5, 5.5)  # Plot margins (t, r, b, l)  (top slightly larger)

############################################################
## Marker display labels (publication)
############################################################
marker_label_map <- c(
  "IL_10"    = "IL-10",
  "IL_1RA"   = "IL-1RA",
  "IL_8_1"   = "IL-8",
  "CCL4"     = "CCL4",
  "CCL3"     = "CCL3",
  "CCL2"     = "CCL2",
  "IL_6"     = "IL-6",
  "IL_1beta" = "IL-1\u03B2",
  "TNF"      = "TNF"
)

############################################################
## 2) Data input
############################################################
studydata <- read_csv("Original_data/studydata_with_ms1_il1r2.csv", show_col_types = FALSE)
wb_long   <- read_csv("Original_data/whole_blood_stimuli_long.csv", show_col_types = FALSE)

############################################################
## 3) Marker order (use your required order)
############################################################
marker_order_raw <- c(
  "IL_10",
  "IL_1RA",
  "IL_8_1",
  "CCL4",
  "CCL3",
  "CCL2",
  "IL_6",
  "IL_1beta",
  "TNF"
)

############################################################
## 4) Core function: Hedges' g + BH (panel-wise)
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
      marker = factor(marker, levels = marker_order_raw)
    ) %>%
    select(marker, value_num, axis_group)
  
  ## Effect size (Hedges' g)  (version-robust)
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
  
  ## t-test + BH within this panel
  tt <- wb_clean %>%
    group_by(marker) %>%
    rstatix::t_test(value_num ~ axis_group, var.equal = FALSE) %>%
    as.data.frame() %>%
    ungroup() %>%
    mutate(p_BH = p.adjust(p, method = "BH")) %>%
    select(marker, statistic, p_BH)
  
  es %>%
    left_join(tt, by = "marker") %>%
    mutate(
      axis = axis_label,
      p_BH.signif = case_when(
        is.na(p_BH) ~ "ns",
        p_BH < 1e-4 ~ "****",
        p_BH < 1e-3 ~ "***",
        p_BH < 1e-2 ~ "**",
        p_BH < 5e-2 ~ "*",
        TRUE ~ "ns"
      ),
      dir = case_when(
        statistic > 0 ~ "UP",
        statistic < 0 ~ "DOWN",
        TRUE ~ NA_character_
      ),
      signif_dir = ifelse(p_BH.signif == "ns" | is.na(dir), "ns", paste(p_BH.signif, dir)),
      marker = factor(as.character(marker), levels = marker_order_raw)
    )
}

############################################################
## 5) Run forest analyses (MS1 and IL1R2)
############################################################
eff_all <- bind_rows(
  compute_hedges_lps(
    studydata, wb_long,
    "MS1_group",
    c("MS1_high", "MS1_low"),
    "MS1: high vs low"
  ),
  compute_hedges_lps(
    studydata, wb_long,
    "IL1R2_group",
    c("IL1R2_positive", "IL1R2_zero"),
    "IL1R2: positive vs zero"
  )
) %>%
  mutate(
    axis   = factor(axis, levels = c("MS1: high vs low", "IL1R2: positive vs zero")),
    marker = factor(marker, levels = marker_order_raw)
  )

############################################################
## 6) Forest plot helpers (zebra, colors, stars)
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
  "* UP"="#ffd4c4","** UP"="#FF9F94","*** UP"="#ff5c33","**** UP"="#cc2900",
  "* DOWN"="#c2d1ff","** DOWN"="#99b3ff","*** DOWN"="#668cff","**** DOWN"="#335CFF"
)

star_df <- eff_all %>%
  filter(p_BH.signif != "ns") %>%
  group_by(axis) %>%
  mutate(
    x_star = max(conf.high, na.rm = TRUE) +
      0.06 * diff(range(conf.low, conf.high, na.rm = TRUE))
  ) %>%
  ungroup()

############################################################
## 6.1) X ranges per facet + top double-headed arrow + labels
############################################################
x_rng <- eff_all %>%
  group_by(axis) %>%
  summarise(
    x_min = min(conf.low,  na.rm = TRUE),
    x_max = max(conf.high, na.rm = TRUE),
    x_pad = 0.12 * (x_max - x_min),
    .groups = "drop"
  ) %>%
  mutate(
    x_left  = x_min - x_pad,
    x_right = x_max + x_pad
  )

arrow_bar_df <- x_rng %>%
  transmute(
    axis,
    xstart = x_left,
    xend   = x_right,
    y      = Inf
  )

arrow_text_df <- x_rng %>%
  mutate(
    x_left_mid  = (x_left + 0) / 2,
    x_right_mid = (0 + x_right) / 2
  ) %>%
  transmute(
    axis,
    x_left_mid, x_right_mid,
    y = Inf,
    lab_left = case_when(
      axis == "MS1: high vs low"        ~ "Lower in MS1 high group",
      axis == "IL1R2: positive vs zero" ~ "Lower in IL1R2 positive group",
      TRUE                              ~ "Lower"
    ),
    lab_right = case_when(
      axis == "MS1: high vs low"        ~ "Higher in MS1 high group",
      axis == "IL1R2: positive vs zero" ~ "Higher in IL1R2 positive group",
      TRUE                              ~ "Higher"
    )
  )

############################################################
## 7) Forest plot (two panels)
############################################################
p_forest <- ggplot(eff_all, aes(x = g, y = marker)) +
  geom_rect(
    data = strip_bg,
    inherit.aes = FALSE,
    aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
    fill = "gray95"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  
  ## Top double-headed arrow (per facet)
  geom_segment(
    data = arrow_bar_df,
    inherit.aes = FALSE,
    aes(x = xstart, xend = xend, y = y, yend = y),
    linewidth = 0.6,
    arrow = arrow(ends = "both", type = "closed", length = unit(3.0, "mm"))
  ) +
  geom_text(
    data = arrow_text_df,
    inherit.aes = FALSE,
    aes(x = x_left_mid, y = y, label = lab_left),
    hjust = 0.5, vjust = 1.25,
    size = (AXIS_TEXT / ggplot2::.pt),
    family = FONT_FAMILY,
    fontface = "bold"
  ) +
  geom_text(
    data = arrow_text_df,
    inherit.aes = FALSE,
    aes(x = x_right_mid, y = y, label = lab_right),
    hjust = 0.5, vjust = 1.25,
    size = (AXIS_TEXT / ggplot2::.pt),
    family = FONT_FAMILY,
    fontface = "bold"
  ) +
  
  geom_segment(aes(x = conf.low, xend = conf.high, yend = marker), linewidth = 0.4) +
  geom_point(aes(fill = signif_dir), shape = 22, size = 4.2, color = "black") +
  geom_text(
    data = star_df,
    inherit.aes = FALSE,
    aes(x = x_star, y = marker, label = p_BH.signif),
    hjust = 0, size = 5,
    family = FONT_FAMILY, color = "black"
  ) +
  scale_fill_manual(values = color_mapping, na.value = "#888888") +
  scale_y_discrete(labels = marker_label_map) +
  facet_grid(cols = vars(axis), scales = "free_x") +
  labs(x = "Hedges’ g (95% CI)", y = NULL) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = BASE_SIZE) +
  theme(
    text = element_text(family = FONT_FAMILY),
    axis.text.y = element_text(face = "bold", size = AXIS_TEXT),
    axis.text.x = element_text(size = AXIS_TEXT),
    axis.title.x = element_text(size = AXIS_TITLE),
    plot.title = element_text(size = TITLE_SIZE),
    
    axis.line.x = element_line(linewidth = LINE_AXIS),
    axis.ticks  = element_line(linewidth = LINE_TICKS),
    
    strip.background = element_blank(),
    strip.text = element_blank(),
    axis.line.y = element_blank(),
    legend.position = "none",
    panel.spacing.x = unit(1.2, "lines"),
    plot.margin = margin(
      t = MARGINS_MM[1], r = MARGINS_MM[2],
      b = MARGINS_MM[3], l = MARGINS_MM[4],
      unit = "mm"
    )
  )

print(p_forest)

ggsave("Figure/LPS_MS1_IL1R2_hedges_final.svg", p_forest, width = 9, height = 8, dpi = 300)


############################################################
## 8) IL1R2 effect WITHIN MS1_high (LPS stimulation)
##    Extension of previous Hedges' g analysis
############################################################

## Subset studydata to MS1_high only
studydata_ms1high <- studydata %>%
  filter(MS1_group == "MS1_high")

## Compute Hedges' g for IL1R2 within MS1_high
eff_ms1high_il1r2 <- compute_hedges_lps(
  studydata = studydata_ms1high,
  wb_long   = wb_long,
  group_var = "IL1R2_group",
  group_levels = c("IL1R2_positive", "IL1R2_zero"),
  axis_label = "IL1R2: positive vs zero (within MS1-high)"
) %>%
  mutate(
    axis   = factor(axis, levels = "IL1R2: positive vs zero (within MS1-high)"),
    marker = factor(marker, levels = marker_order_raw)
  )

## Zebra background
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

## Star positions (consistent logic)
star_df_ms1high <- eff_ms1high_il1r2 %>%
  filter(p_BH.signif != "ns") %>%
  group_by(axis) %>%
  mutate(
    x_star = max(conf.high, na.rm = TRUE) +
      0.06 * diff(range(conf.low, conf.high, na.rm = TRUE))
  ) %>%
  ungroup()

## X range + top arrow + labels (single panel)
x_rng_ms1high <- eff_ms1high_il1r2 %>%
  summarise(
    x_min = min(conf.low,  na.rm = TRUE),
    x_max = max(conf.high, na.rm = TRUE)
  ) %>%
  mutate(
    x_pad   = 0.12 * (x_max - x_min),
    x_left  = x_min - x_pad,
    x_right = x_max + x_pad
  )

arrow_bar_ms1high <- x_rng_ms1high %>%
  transmute(
    xstart = x_left,
    xend   = x_right,
    y      = Inf
  )

arrow_text_ms1high <- x_rng_ms1high %>%
  transmute(
    x_left_mid  = (x_left + 0) / 2,
    x_right_mid = (0 + x_right) / 2,
    y = Inf,
    lab_left  = "Lower in IL1R2 positive group",
    lab_right = "Higher in IL1R2 positive group"
  )

## Plot: IL1R2 within MS1_high (match GLOBAL style)
p_ms1high_il1r2 <- ggplot(eff_ms1high_il1r2, aes(x = g, y = marker)) +
  geom_rect(
    data = strip_bg_ms1high,
    inherit.aes = FALSE,
    aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
    fill = "gray95"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  
  ## Top double-headed arrow + bold labels (like your reference)
  geom_segment(
    data = arrow_bar_ms1high,
    inherit.aes = FALSE,
    aes(x = xstart, xend = xend, y = y, yend = y),
    linewidth = 0.6,
    arrow = arrow(ends = "both", type = "closed", length = unit(3.0, "mm"))
  ) +
  geom_text(
    data = arrow_text_ms1high,
    inherit.aes = FALSE,
    aes(x = x_left_mid, y = y, label = lab_left),
    hjust = 0.5, vjust = 1.25,
    size = 4.5,
    family = FONT_FAMILY,
    fontface = "bold"
  ) +
  geom_text(
    data = arrow_text_ms1high,
    inherit.aes = FALSE,
    aes(x = x_right_mid, y = y, label = lab_right),
    hjust = 0.5, vjust = 1.25,
    size = 4.5,
    family = FONT_FAMILY,
    fontface = "bold"
  ) +
  
  geom_segment(aes(x = conf.low, xend = conf.high, yend = marker), linewidth = 0.4) +
  geom_point(aes(fill = signif_dir), shape = 22, size = 4.2, color = "black") +
  
  ## Stars (same as your fixed code)
  geom_text(
    data = star_df_ms1high,
    inherit.aes = FALSE,
    aes(x = x_star, y = marker, label = p_BH.signif),
    hjust = 0, size = 5,
    family = FONT_FAMILY, color = "black"
  ) +
  
  scale_fill_manual(values = color_mapping, na.value = "#888888") +
  scale_y_discrete(labels = marker_label_map) +
  facet_grid(cols = vars(axis), scales = "free_x") +
  labs(x = "Hedges’ g (95% CI)", y = NULL) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = BASE_SIZE) +
  theme(
    text = element_text(family = FONT_FAMILY),
    axis.text.y = element_text(face = "bold", size = AXIS_TEXT),
    axis.text.x = element_text(size = AXIS_TEXT),
    axis.title.x = element_text(size = AXIS_TITLE),
    plot.title = element_text(size = TITLE_SIZE),
    
    axis.line.x = element_line(linewidth = LINE_AXIS),
    axis.ticks  = element_line(linewidth = LINE_TICKS),
    
    strip.background = element_blank(),
    strip.text = element_blank(),
    axis.line.y = element_blank(),
    legend.position = "none",
    panel.spacing.x = unit(1.2, "lines"),
    plot.margin = margin(
      t = MARGINS_MM[1], r = MARGINS_MM[2],
      b = MARGINS_MM[3], l = MARGINS_MM[4],
      unit = "mm"
    )
  )

print(p_ms1high_il1r2)

ggsave(
  "Figure/LPS_IL1R2_within_MS1_high_hedges.svg",
  p_ms1high_il1r2,
  width = 5,
  height = 7,
  dpi = 300
)


############################################################
## 9) TNF: Panel A (MS1 spline scatter) + Panel B (boxplot + curve)
##    Requirement:
##    - B keeps boxplot (0..6 groups)
##    - draw a smooth curve ON TOP of the boxplot
##    - curve uses SAME method as A (GAM)
##    - curve starts from group 1 (exclude group 0)
############################################################

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
MARGINS_MM  <- c(5.5, 5.5, 5.5, 5.5)  # Plot margins (t, r, b, l)

############################################################
## Theme (GLOBAL style)
############################################################
base_theme <- theme_classic(base_size = BASE_SIZE) +
  theme(
    text        = element_text(family = FONT_FAMILY),
    axis.title  = element_text(size = AXIS_TITLE),
    axis.text   = element_text(size = AXIS_TEXT),
    plot.title  = element_text(size = TITLE_SIZE),
    axis.line   = element_line(linewidth = LINE_AXIS),
    axis.ticks  = element_line(linewidth = LINE_TICKS),
    plot.margin = margin(
      t = MARGINS_MM[1], r = MARGINS_MM[2],
      b = MARGINS_MM[3], l = MARGINS_MM[4],
      unit = "mm"
    )
  )

col_ms1   <- "#2F4B7C"
col_il1r2 <- "#A23B3B"

############################################################
## Prepare TNF (LPS)
############################################################
df_tnf <- wb_long %>%
  filter(stimulation == "LPS", marker == "TNF") %>%
  mutate(value = as.numeric(value)) %>%
  left_join(
    studydata %>% select(EB_id, MS1_percent, IL1R2_percent),
    by = c("ID" = "EB_id")
  ) %>%
  filter(!is.na(value), !is.na(MS1_percent), !is.na(IL1R2_percent)) %>%
  mutate(
    value_log = log10(pmax(value, 1e-9)),
    IL1R2 = IL1R2_percent
  )

############################################################
## Panel A (MS1 smooth)
############################################################
p_ms1 <- ggplot(df_tnf, aes(x = MS1_percent, y = value_log)) +
  geom_point(color = col_ms1, alpha = 0.6, size = 2) +
  geom_smooth(
    method  = "gam",
    formula = y ~ s(x, k = 4),
    color   = col_ms1,
    fill    = col_ms1,
    alpha   = 0.25,
    linewidth = 1
  ) +
  labs(
    x = "Percentage of MS1 cells",
    y = "LPS-induced TNF (log10)"
  ) +
  base_theme

############################################################
## Build IL1R2 groups: 0 (==0) + 1..6 among >0
############################################################
df_il1r2_grp <- df_tnf %>%
  mutate(IL1R2_pos = IL1R2 > 0)

pos_breaks <- quantile(
  df_il1r2_grp$IL1R2[df_il1r2_grp$IL1R2 > 0],
  probs = seq(0, 1, length.out = 7),
  na.rm = TRUE
)

pos_breaks <- unique(as.numeric(pos_breaks))
if (length(pos_breaks) < 7) {
  pos_breaks <- pretty(
    range(df_il1r2_grp$IL1R2[df_il1r2_grp$IL1R2 > 0], na.rm = TRUE),
    n = 6
  )
}

df_il1r2_grp <- df_il1r2_grp %>%
  mutate(
    IL1R2_group7 = case_when(
      IL1R2 == 0 ~ 0L,
      IL1R2 > 0  ~ as.integer(cut(
        IL1R2,
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
## Fit GAM curve ONLY from group 1..6 (start from 1)
############################################################
df_for_curve <- df_il1r2_grp %>% filter(IL1R2_group7_num >= 1)

gam_grp <- mgcv::gam(value_log ~ s(IL1R2_group7_num, k = 4), data = df_for_curve)

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
## Panel B: boxplot + jitter + smooth curve (same GAM style)
############################################################
p_il1r2_box_curve <- ggplot(df_il1r2_grp, aes(x = IL1R2_group7_num, y = value_log)) +
  geom_boxplot(
    aes(group = IL1R2_group7_num),
    outlier.shape = NA,
    linewidth = 0.8,
    fill = "#F3D7CF"
  ) +
  geom_jitter(
    width = 0.18,
    height = 0,
    alpha = 0.55,
    size = 2,
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
    linewidth = 1
  ) +
  scale_x_continuous(
    breaks = 0:6,
    labels = 0:6
  ) +
  labs(
    x = "IL1R2+ cells group",
    y = "LPS-induced TNF (log10)"
  ) +
  base_theme

############################################################
## Combine A + B with A/B tags
############################################################
p_tnf_final <- p_ms1 + p_il1r2_box_curve +
  patchwork::plot_layout(ncol = 2) +
  patchwork::plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(
        family = FONT_FAMILY,
        size = TAG_SIZE,
        face = "plain"
      )
    )
  )

print(p_tnf_final)

ggsave(
  "Figure/TNF_MS1_A_IL1R2_B_boxplot_with_curve.svg",
  p_tnf_final,
  width = 10,
  height = 4.5,
  dpi = 300
)