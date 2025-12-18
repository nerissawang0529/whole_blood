rm(list = ls())

## install.packages("pacman")
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
               tableone, kmed, survival, survminer, skimr, lattice, reshape2, grDevices, ggrepel)

studydata_with_cts <- read.csv("Original_data/studydata_with_cts.csv", check.names = FALSE)
whole_blood_stimuli_long <- read.csv("Original_data/whole_blood_stimuli_long.csv", check.names = FALSE)
whole_blood_stimuli_long$HiIL6 <- NULL
whole_blood_stimuli_long$level <- NULL

write.csv(whole_blood_stimuli_long, "Original_data/stimulation_EB.csv", row.names = FALSE)








# Make sure IDs are comparable
id_key <- studydata_with_cts %>%
  transmute(EB_id = as.character(EB_id),
            CTS   = as.character(CTS))

wb <- whole_blood_stimuli_long %>%
  mutate(ID = as.character(ID)) %>%
  semi_join(id_key, by = c("ID" = "EB_id")) %>%      # keep EB_id only
  left_join(id_key, by = c("ID" = "EB_id"))          # add CTS column

## ===== Packages =====
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr); library(forcats); library(rlang)
  library(rstatix); library(purrr); library(tidyr); library(ggplot2)
})

## ===== Settings =====
stim_keep <- c("KP","PP","LPS")
group_var <- "CTS"   # groups: 1,2,3
analysis_scale <- "log10"  # or "raw"

## ===== 1) Clean & subset =====
wb_clean <- wb %>%
  dplyr::filter(!is.na(marker), !is.na(stimulation), !is.na(.data[[group_var]])) %>%
  dplyr::filter(stimulation %in% stim_keep) %>%
  dplyr::mutate(
    value_raw = as.character(value),
    v_dot   = readr::parse_number(value_raw, locale = readr::locale(decimal_mark = ".", grouping_mark = ",")),
    v_comma = readr::parse_number(value_raw, locale = readr::locale(decimal_mark = ",", grouping_mark = ".")),
    value   = dplyr::coalesce(v_dot, v_comma),
    value   = dplyr::if_else(stringr::str_detect(value_raw, "^\\s*<\\s*\\d+(?:[\\.,]\\d+)?\\s*$"),
                             readr::parse_number(stringr::str_replace(value_raw, "<","")), value),
    value   = dplyr::if_else(stringr::str_detect(value_raw, "^\\s*>\\s*\\d+(?:[\\.,]\\d+)?\\s*$"),
                             readr::parse_number(stringr::str_replace(value_raw, ">","")), value),
    value_num = if (analysis_scale == "log10") log10(pmax(value, 1e-9)) else value,
    stimulation = factor(stimulation, levels = stim_keep),
    !!group_var := factor(.data[[group_var]])
  ) %>%
  dplyr::select(stimulation, marker, value_num, !!sym(group_var))

## ===== 2) Specified contrasts (positive = higher in first group) =====
pairs_df <- tibble::tribble(
  ~first, ~second,
  "1",    "2",
  "2",    "3",
  "1",    "3"
)

effsize_df <- purrr::pmap_dfr(pairs_df, function(first, second) {
  dat <- wb_clean %>%
    dplyr::filter(.data[[group_var]] %in% c(first, second)) %>%
    dplyr::mutate(!!group_var := forcats::fct_relevel(.data[[group_var]], second, first))
  
  es <- dat %>%
    dplyr::group_by(stimulation, marker) %>%
    rstatix::cohens_d(reformulate(group_var, response = "value_num"),
                      var.equal = FALSE, hedges.correction = TRUE, ci = TRUE) %>%
    dplyr::ungroup() %>%
    dplyr::rename(g = effsize)
  
  tt <- dat %>%
    dplyr::group_by(stimulation, marker) %>%
    rstatix::t_test(reformulate(group_var, response = "value_num"), var.equal = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(p_BH = p.adjust(p, method = "BH")) %>%
    dplyr::select(stimulation, marker, statistic, p_BH)
  
  es %>%
    dplyr::left_join(tt, by = c("stimulation","marker")) %>%
    dplyr::mutate(
      contrast_lab = sprintf("CTS%s vs CTS%s\n(positive = higher in first group)", first, second),
      p_BH.signif = dplyr::case_when(
        is.na(p_BH) ~ "ns",
        p_BH < 1e-4 ~ "****",
        p_BH < 1e-3 ~ "***",
        p_BH < 1e-2 ~ "**",
        p_BH < 5e-2 ~ "*",
        TRUE ~ "ns"
      ),
      dir = dplyr::case_when(statistic > 0 ~ "UP",
                             statistic < 0 ~ "DOWN",
                             TRUE ~ NA_character_),
      signif_dir = ifelse(p_BH.signif == "ns" | is.na(dir), "ns", paste(p_BH.signif, dir))
    )
})

## ===== 3) Pretty labels =====
stim_labels <- c(KP = "Kleb.", PP = "Pneumo.", LPS = "LPS")
marker_pretty <- c(
  "IL_1beta"="IL-1β","IL_1RA"="IL-1RA","IL_8_1"="IL-8","IL_6"="IL-6",
  "IL_10"="IL-10","TNF"="TNF","CCL2"="CCL2","CCL3"="CCL3","CCL4"="CCL4"
)

effsize_df <- effsize_df %>%
  dplyr::mutate(
    stimulation = forcats::fct_relabel(stimulation, ~ stim_labels[.x] %||% .x),
    marker = factor(marker_pretty[as.character(marker)] %||% as.character(marker))
  )

marker_order <- c("CCL2","CCL3","CCL4","IL-1β","IL-6","TNF","IL-8","IL-1RA","IL-10")
contrast_order <- c(
  "CTS1 vs CTS2\n(positive = higher in first group)",
  "CTS2 vs CTS3\n(positive = higher in first group)",
  "CTS1 vs CTS3\n(positive = higher in first group)"
)

effsize_df <- effsize_df %>%
  dplyr::mutate(
    marker = factor(marker, levels = marker_order),
    contrast_lab = factor(contrast_lab, levels = contrast_order)
  )

## ===== 4) Striped backgrounds =====
strip_bg <- effsize_df %>%
  dplyr::group_by(stimulation, contrast_lab) %>%
  dplyr::distinct(marker) %>%
  dplyr::arrange(stimulation, contrast_lab, marker) %>%
  dplyr::mutate(row_id = dplyr::row_number(),
                ymin = row_id - 0.5, ymax = row_id + 0.5) %>%
  dplyr::filter(row_id %% 2 == 1) %>%
  dplyr::ungroup()

## ===== 5) Label map for new column titles =====
lab_map <- c(
  "CTS1 vs CTS2\n(positive = higher in first group)" = 
    "CTS2 vs CTS1\n(positive = higher in first group)",
  "CTS2 vs CTS3\n(positive = higher in first group)" = 
    "CTS3 vs CTS2\n(positive = higher in first group)",
  "CTS1 vs CTS3\n(positive = higher in first group)" = 
    "CTS3 vs CTS1\n(positive = higher in first group)"
)

## ===== 6) Plot =====
color_mapping <- c(
  "ns"="#888888",
  "* UP"="#ffd4c4","** UP"="#FF9F94","*** UP"="#ff5c33","**** UP"="#cc2900",
  "* DOWN"="#c2d1ff","** DOWN"="#99b3ff","*** DOWN"="#668cff","**** DOWN"="blue"
)

library(ggplot2)
library(dplyr)
library(grid)

p_all <- ggplot(effsize_df, aes(x = g, y = marker)) +
  geom_rect(
    data = strip_bg, inherit.aes = FALSE,
    aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
    fill = "gray95", colour = NA
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.4) +
  geom_segment(aes(x = conf.low, xend = conf.high, yend = marker), size = 0.4) +
  geom_point(aes(fill = signif_dir), shape = 22, size = 3) +
  geom_text(
    data = effsize_df %>%
      group_by(stimulation, contrast_lab) %>%
      mutate(x_lab = max(conf.high, na.rm = TRUE) + 0.10) %>%
      ungroup(),
    aes(x = x_lab, label = ifelse(p_BH.signif == "ns", "", p_BH.signif)),
    hjust = 0, size = 3.3, family = "sans"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.35))) +
  scale_y_discrete(drop = TRUE, expand = expansion(add = c(0.5, 0.5))) +
  scale_fill_manual(values = color_mapping, na.value = "#888888") +
  facet_grid(
    rows = vars(stimulation),
    cols = vars(contrast_lab),
    scales = "free_x",
    labeller = labeller(contrast_lab = as_labeller(lab_map))
  ) +
  labs(
    x = "Hedges’ g (positive = higher in first group), 95% CI",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y     = element_text(face = "bold"),
    strip.text.x    = element_text(face = "bold"),
    strip.text.y    = element_text(face = "bold"),
    panel.grid      = element_blank(),
    axis.line.x     = element_line(size = 0.25),
    legend.position = "none",
    plot.margin     = ggplot2::margin(10, 30, 10, 10, unit = "pt"),
    panel.spacing.x = grid::unit(10, "pt"),
    panel.spacing.y = grid::unit(8, "pt"),
    plot.caption    = element_text(hjust = 0, size = 10, color = "gray30")
  )


print(p_all)

ggsave("Figure/cts_hedges_stimuli.svg",
       p_all, width = 12, height = 8, dpi = 300)


## ===== 依赖包 =====
library(dplyr)
library(ggplot2)
library(forcats)
library(patchwork)

## ===== 统一配色 =====
colors <- c("1"="#f28e8e",  # CTS1：珊瑚红
            "2"="#5bc0de",  # CTS2：天蓝
            "3"="#a1d99b")  # CTS3：草绿

## ===== 通用绘图函数 =====
make_plot <- function(stim_now, first, second, contrast_lab, stim_label) {
  dat_sub <- wb_clean %>%
    filter(stimulation == stim_now, CTS %in% c(first, second)) %>%
    mutate(CTS = fct_relevel(CTS, first, second))
  
  ggplot(dat_sub, aes(x = CTS, y = value_num, fill = CTS)) +
    geom_boxplot(width = 0.6, outlier.size = 0.8) +
    geom_jitter(width = 0.1, alpha = 0.35, size = 0.8) +
    facet_wrap(~ marker, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = colors, drop = FALSE) +
    labs(
      title = paste0(stim_label, " stimulation — ", contrast_lab),
      x = NULL,
      y = "log10(value)"
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}

## ======== Kleb. (KP) ========
p1_kp <- make_plot("KP", "2", "1", "CTS2 vs CTS1", "Kleb.")
p2_kp <- make_plot("KP", "3", "2", "CTS3 vs CTS2", "Kleb.")
p3_kp <- make_plot("KP", "3", "1", "CTS3 vs CTS1", "Kleb.")

final_plot_kp <- p1_kp / p2_kp / p3_kp +
  plot_annotation(
    title = "Kleb. stimulation — CTS contrasts across all markers",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  )

ggsave("Figure/Box_allmarkers_Kleb_stimuli.svg", final_plot_kp, width = 12, height = 8, dpi = 300)
print(final_plot_kp)


## ======== Pneumo. (PP) ========
p1_pp <- make_plot("PP", "2", "1", "CTS2 vs CTS1", "Pneumo.")
p2_pp <- make_plot("PP", "3", "2", "CTS3 vs CTS2", "Pneumo.")
p3_pp <- make_plot("PP", "3", "1", "CTS3 vs CTS1", "Pneumo.")

final_plot_pp <- p1_pp / p2_pp / p3_pp +
  plot_annotation(
    title = "Pneumo. stimulation — CTS contrasts across all markers",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  )

ggsave("Figure/Box_allmarkers_Pneumo_stimuli.svg", final_plot_pp, width = 12, height = 8, dpi = 300)
print(final_plot_pp)


## ======== LPS ========
p1_lps <- make_plot("LPS", "2", "1", "CTS2 vs CTS1", "LPS")
p2_lps <- make_plot("LPS", "3", "2", "CTS3 vs CTS2", "LPS")
p3_lps <- make_plot("LPS", "3", "1", "CTS3 vs CTS1", "LPS")

final_plot_lps <- p1_lps / p2_lps / p3_lps +
  plot_annotation(
    title = "LPS stimulation — CTS contrasts across all markers",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  )

ggsave("Figure/Box_allmarkers_LPS_stimuli.svg", final_plot_lps, width = 12, height = 8, dpi = 300)
print(final_plot_lps)
