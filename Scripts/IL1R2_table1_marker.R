rm(list = ls ())

## =========================
## 0) Load dependencies (must be installed already)
## =========================
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
               tableone, kmed, survival, survminer, skimr, lattice, reshape2, grDevices, ggrepel)



## ===============================
## IL1R2 grouping → Table 1 + Hedges’ g
## ===============================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2); library(tableone); library(forcats)
  library(tidyr); library(purrr); library(rstatix); library(ggpubr)
})

# 1️⃣ 读取 CIBERSORTx 输出并选择 IL1R2 列
il1r2_df <- readr::read_csv("Original_data/CIBERSORTx_IL1R2_percent.csv",
                            show_col_types = FALSE)

il1r2_df <- il1r2_df %>%
  dplyr::rename(Sample = Mixture, IL1R2 = `IL1R2_percent`) %>%
  dplyr::mutate(IL1R2 = as.numeric(IL1R2))

stopifnot(!anyNA(il1r2_df$IL1R2))

# 2️⃣ Map Sample → ID
dfwb_sel <- readr::read_csv("Original_data/dfwb_sel.csv", show_col_types = FALSE)
il1r2_map <- il1r2_df %>%
  left_join(dfwb_sel %>% select(Sample, ID), by = "Sample") %>%
  filter(!is.na(ID))

# 3️⃣ 分组逻辑：0 / low / high（非零样本按上四分位分）
p75_nonzero <- quantile(il1r2_map$IL1R2[il1r2_map$IL1R2 > 0], 0.75, na.rm = TRUE)

il1r2_map <- il1r2_map %>%
  mutate(
    IL1R2_group = case_when(
      IL1R2 == 0 ~ "IL1R2⁻",
      IL1R2 <= p75_nonzero ~ "IL1R2⁺ low",
      IL1R2 > p75_nonzero ~ "IL1R2⁺ high"
    )
  ) %>%
  mutate(IL1R2_group = factor(IL1R2_group, levels = c("IL1R2⁻","IL1R2⁺ low","IL1R2⁺ high")))

table(il1r2_map$IL1R2_group)

# 4️⃣ 合并临床数据
studydata <- readr::read_csv("Original_data/studydata.csv", show_col_types = FALSE)

studydata_with_il1r2 <- studydata %>%
  mutate(EB_id = as.character(EB_id)) %>%
  inner_join(il1r2_map %>% transmute(EB_id = as.character(ID),
                                     Sample,
                                     IL1R2, IL1R2_group),
             by = "EB_id")

# 保存临床合并文件
write.csv(studydata_with_il1r2,
          "Original_data/studydata_with_il1r2.csv",
          row.names = FALSE)
cat("✅ Saved: Original_data/studydata_with_il1r2.csv\n")


## ===============================
##  Table 1
## ===============================
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
             "hospdeath",	"mortality_d30", "mortality_d90","IL1R2_group")

catvars <- c("group","flow_Data", "spectral", "inclusion_hospital", "gender",	"ethnic_group#White",	
             "hypertension",	"cpd",	"COPD",	"diabetes",	
             "ccd","ckd",	"mneoplasm",	"immune_sup",	"cnd",	
             "Oxygen_therapy_1", 
             "antibiotic_seven_days", 
             "ICU_Medium_Care_admission_1",	 
             "Non_invasive_ventilation_1",	"Invasive_ventilation_1",	
             "bacteremia",
             "hospdeath",	"mortality_d30", "mortality_d90", "pathogen_cultured","IL1R2_group")

nonnormal <- c("sampling_time",		"age_yrs", "BMI", "symptom_days",  "CCI", "qSOFA_score",	"MEWS_score",	"CURB_score",	"PSI_new",
               "crp_1_1",	 "WBC_2_1",	"Neutrophil_unit_1",	"Lymphocyte_1_1",	
               "Creatinine_value_1",
               "Platelets_value_1",	
               "Blood_Urea_Nitrogen_value_1",	"length_of_oxygen","length_of_stay" ,"icu_stay", "lenght_of_intubation")


tab1 <- CreateTableOne(vars = allvars, data = studydata_with_il1r2,
                       strata = "IL1R2_group", factorVars = catvars, test = TRUE)
tab1_df <- as.data.frame(print(tab1, nonnormal = nonnormal, quote = FALSE, noSpaces = TRUE, smd = TRUE, missing = TRUE))
tab1_df$tablename <- rownames(tab1_df)
tab1_df <- tab1_df[, c("tablename", setdiff(names(tab1_df), "tablename"))]
write.csv(tab1_df, "Original_data/table1_il1r2.csv", row.names = FALSE)
cat("✅ Table 1 exported to: Original_data/table1_il1r2.csv\n")



## ===============================
##  Hedges’ g analysis
## ===============================
whole_blood_stimuli_long <- read_csv("Original_data/whole_blood_stimuli_long.csv", show_col_types = FALSE)
if ("HiIL6" %in% names(whole_blood_stimuli_long)) whole_blood_stimuli_long$HiIL6 <- NULL

id_key <- studydata_with_il1r2 %>%
  transmute(EB_id = as.character(EB_id),
            Sample = as.character(Sample),
            IL1R2 = as.numeric(IL1R2),
            IL1R2_group = factor(IL1R2_group,
                                 levels = c("IL1R2⁻","IL1R2⁺ low","IL1R2⁺ high")))

wb <- whole_blood_stimuli_long %>%
  mutate(ID = as.character(ID)) %>%
  semi_join(id_key, by = c("ID" = "EB_id")) %>%
  left_join(id_key, by = c("ID" = "EB_id"))

stim_keep <- c("KP","PP","LPS")
group_var <- "IL1R2_group"
analysis_scale <- "log10"

wb_clean <- wb %>%
  filter(!is.na(marker), !is.na(stimulation), !is.na(.data[[group_var]])) %>%
  filter(stimulation %in% stim_keep) %>%
  mutate(
    value_raw = as.character(value),
    value     = readr::parse_number(value_raw),
    value_num = if (analysis_scale == "log10") log10(pmax(value, 1e-9)) else value,
    stimulation = factor(stimulation, levels = stim_keep),
    !!group_var := factor(.data[[group_var]], levels = c("IL1R2⁻","IL1R2⁺ low","IL1R2⁺ high"))
  ) %>%
  select(stimulation, marker, value_num, !!sym(group_var))

# === Effect size 计算保持不变 ===
il1_lvls <- levels(wb_clean[[group_var]])
contrast_df <- as.data.frame(t(combn(il1_lvls, 2)))
names(contrast_df) <- c("ref","comp")

effsize_df <- purrr::pmap_dfr(contrast_df, function(ref, comp) {
  dat <- wb_clean %>%
    filter(.data[[group_var]] %in% c(ref, comp)) %>%
    mutate(!!group_var := forcats::fct_relevel(.data[[group_var]], ref, comp))
  
  es <- dat %>%
    group_by(stimulation, marker) %>%
    rstatix::cohens_d(reformulate(group_var, response = "value_num"),
                      var.equal = FALSE, hedges.correction = TRUE, ci = TRUE) %>%
    ungroup() %>%
    rename_with(~"g", matches("effsize"))   
  
  tt <- dat %>%
    group_by(stimulation, marker) %>%
    rstatix::t_test(reformulate(group_var, response = "value_num"), var.equal = FALSE) %>%
    ungroup() %>%
    mutate(p_BH = p.adjust(p, method = "BH")) %>%
    select(stimulation, marker, statistic, p_BH)
  
  es %>%
    left_join(tt, by = c("stimulation","marker")) %>%
    mutate(
      contrast_lab = sprintf("%s vs %s (ref = %s)", comp, ref, ref),
      p_BH.signif = case_when(
        is.na(p_BH) ~ "ns",
        p_BH < 1e-4 ~ "****",
        p_BH < 1e-3 ~ "***",
        p_BH < 1e-2 ~ "**",
        p_BH < 5e-2 ~ "*",
        TRUE ~ "ns"
      ),
      dir = case_when(statistic > 0 ~ "UP", statistic < 0 ~ "DOWN", TRUE ~ NA_character_),
      signif_dir = ifelse(p_BH.signif == "ns" | is.na(dir), "ns", paste(p_BH.signif, dir))
    )
})

# === Facet label 显示修改（仅改标题）===
lab_map <- c(
  "IL1R2⁺ high vs IL1R2⁻ (ref = IL1R2⁻)" = "IL1R2⁻ vs IL1R2⁺ high (positive = higher in first group)",
  "IL1R2⁺ low vs IL1R2⁻ (ref = IL1R2⁻)"  = "IL1R2⁻ vs IL1R2⁺ low (positive = higher in first group)",
  "IL1R2⁺ high vs IL1R2⁺ low (ref = IL1R2⁺ low)" = "IL1R2⁺ low vs IL1R2⁺ high (positive = higher in first group)"
)

stim_labels <- c(KP = "Kleb.", PP = "Pneumo.", LPS = "LPS")
marker_pretty <- c("IL_1beta"="IL-1β","IL_1RA"="IL-1RA","IL_8_1"="IL-8","IL_6"="IL-6",
                   "IL_10"="IL-10","TNF"="TNF","CCL2"="CCL2","CCL3"="CCL3","CCL4"="CCL4")

marker_order <- c("CCL2","CCL3","CCL4","IL-1β","IL-6","TNF","IL-8","IL-1RA","IL-10")
contrast_order <- c(
  "IL1R2⁺ high vs IL1R2⁻ (ref = IL1R2⁻)",
  "IL1R2⁺ low vs IL1R2⁻ (ref = IL1R2⁻)",
  "IL1R2⁺ high vs IL1R2⁺ low (ref = IL1R2⁺ low)"
)

effsize_df <- effsize_df %>%
  mutate(
    stimulation = fct_relabel(stimulation, ~ stim_labels[.x] %||% .x),
    marker = factor(marker_pretty[as.character(marker)] %||% as.character(marker),
                    levels = marker_order),
    contrast_lab = factor(contrast_lab, levels = contrast_order)
  )

# === 灰条背景 ===
strip_bg <- effsize_df %>%
  group_by(stimulation, contrast_lab) %>%
  distinct(marker) %>%
  arrange(stimulation, contrast_lab, marker) %>%
  mutate(row_id = row_number(),
         ymin = row_id - 0.5, ymax = row_id + 0.5) %>%
  filter(row_id %% 2 == 1) %>%
  ungroup()

color_mapping <- c(
  "ns"="#888888",
  "* UP"="#ffd4c4","** UP"="#FF9F94","*** UP"="#ff5c33","**** UP"="#cc2900",
  "* DOWN"="#c2d1ff","** DOWN"="#99b3ff","*** DOWN"="#668cff","**** DOWN"="blue"
)

# === 绘图 ===
p_il1r2_all <- ggplot(effsize_df, aes(x = g, y = marker)) +
  geom_rect(data = strip_bg, inherit.aes = FALSE,
            aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
            fill = "gray95", colour = NA) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_segment(aes(x = conf.low, xend = conf.high, yend = marker), linewidth = 0.4) +
  geom_point(aes(fill = signif_dir), shape = 22, size = 3) +
  scale_fill_manual(values = color_mapping, na.value = "#888888") +
  facet_grid(rows = vars(stimulation),
             cols = vars(contrast_lab),
             scales = "free_x",
             labeller = labeller(contrast_lab = as_labeller(lab_map))) +
  labs(x = "Hedges’ g (positive = higher in first group)", y = NULL,
       title = "Stimulated whole blood: IL1R2-group contrasts") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y   = element_text(face = "bold"),
    strip.text.x  = element_text(face = "bold"),
    strip.text.y  = element_text(face = "bold"),
    panel.grid    = element_blank(),
    legend.position = "none",
    plot.margin   = margin(10, 30, 10, 10)
  )

dir.create("Figure", showWarnings = FALSE)
ggsave("Figure/il1r2_hedges_labelReversed.svg", p_il1r2_all, width = 12, height = 8, dpi = 300)
cat("✅ Figure saved: Figure/il1r2_hedges_labelReversed.svg\n")


## ===============================
##  Box plot（与 Hedges’ g 对应的标注）
## ===============================
stim_now <- "KP"
dat_kleb <- wb_clean %>%
  filter(stimulation == stim_now) %>%
  mutate(marker = factor(marker),
         IL1R2_group = as.factor(IL1R2_group))

colors <- c("IL1R2⁻"="#f28e8e", "IL1R2⁺ low"="#5bc0de", "IL1R2⁺ high"="#a1d99b")

make_plot <- function(first, second, contrast_lab) {
  dat_sub <- dat_kleb %>%
    filter(IL1R2_group %in% c(first, second)) %>%
    mutate(IL1R2_group = fct_relevel(IL1R2_group, first, second))
  
  ggplot(dat_sub, aes(x = IL1R2_group, y = value_num, fill = IL1R2_group)) +
    geom_boxplot(width = 0.6, outlier.size = 0.8) +
    geom_jitter(width = 0.1, alpha = 0.35, size = 0.8) +
    facet_wrap(~ marker, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = colors, drop = FALSE) +
    labs(
      title = paste0("Kleb. stimulation — ", contrast_lab),
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

p1 <- make_plot("IL1R2⁻", "IL1R2⁺ high", "IL1R2⁻ vs IL1R2⁺ high")
p2 <- make_plot("IL1R2⁻", "IL1R2⁺ low",  "IL1R2⁻ vs IL1R2⁺ low")
p3 <- make_plot("IL1R2⁺ low", "IL1R2⁺ high", "IL1R2⁺ low vs IL1R2⁺ high")

final_plot <- p1 / p2 / p3 +
  plot_annotation(
    title = "Kleb. stimulation — IL1R2 contrasts across all markers",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  )

print(final_plot)
ggsave("Figure/Box_allmarkers_IL1R2_labelReversed.svg", final_plot, width = 12, height = 8, dpi = 300)
cat("✅ Figure saved: Figure/Box_allmarkers_IL1R2_labelReversed.svg\n")



## ===============================
##  IL1R2 spline plots (lm + rcs)
## ===============================

library(dplyr)
library(ggplot2)
library(rms)
library(readr)
library(patchwork)

# ---- Load merged IL1R2 + clinical data ----
studydata_with_il1r2 <- read_csv("Original_data/studydata_with_il1r2.csv", show_col_types = FALSE)
whole_blood_stimuli_long <- read_csv("Original_data/whole_blood_stimuli_long.csv", show_col_types = FALSE)

# ---- Markers & stimulations of interest ----
markers_keep <- c("TNF", "IL_1RA")
stim_keep <- c("LPS", "PP", "KP")

# ---- Prepare dataset ----
wb <- whole_blood_stimuli_long %>%
  filter(marker %in% markers_keep, stimulation %in% stim_keep) %>%
  mutate(ID = as.character(ID)) %>%
  left_join(
    studydata_with_il1r2 %>%
      transmute(EB_id = as.character(EB_id), IL1R2 = as.numeric(IL1R2)),
    by = c("ID" = "EB_id")
  ) %>%
  mutate(
    value = dplyr::coalesce(
      suppressWarnings(as.numeric(value)),
      readr::parse_number(as.character(value))
    ),
    log10_value = log10(pmax(value, 1e-9))
  ) %>%
  filter(!is.na(IL1R2), !is.na(log10_value))

# ---- Helper function: fit & plot spline ----
plot_spline <- function(stim_now, marker_now) {
  dat <- wb %>% filter(stimulation == stim_now, marker == marker_now)
  
  # linear model with restricted cubic spline (3 knots)
  fit <- lm(log10_value ~ rcs(IL1R2, 3), data = dat)
  
  # predicted line and 95% CI
  newdat <- data.frame(IL1R2 = seq(min(dat$IL1R2), max(dat$IL1R2), length.out = 200))
  p <- predict(fit, newdata = newdat, se.fit = TRUE)
  newdat$fit   <- p$fit
  newdat$upper <- p$fit + 1.96 * p$se.fit
  newdat$lower <- p$fit - 1.96 * p$se.fit
  
  ggplot(newdat, aes(x = IL1R2, y = fit)) +
    geom_line(color = "#0072B2", linewidth = 1) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#0072B2", alpha = 0.2) +
    labs(
      title = paste(stim_now, "—", marker_now),
      x = "IL1R2 score",
      y = paste0("log10(", marker_now, ")")
    ) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

# ---- Make all six plots ----
plots <- list(
  plot_spline("LPS", "TNF"),
  plot_spline("PP",  "TNF"),
  plot_spline("KP",  "TNF"),
  plot_spline("LPS", "IL_1RA"),
  plot_spline("PP",  "IL_1RA"),
  plot_spline("KP",  "IL_1RA")
)

# ---- Arrange vertically: same marker in one column ----
final_plot <- (plots[[1]] / plots[[2]] / plots[[3]]) | 
  (plots[[4]] / plots[[5]] / plots[[6]]) +
  plot_annotation(
    title = "Associations between IL1R2 score and cytokine levels",
    subtitle = "Modeled by lm(log10(marker) ~ rcs(IL1R2, 3))",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  )

# ---- Save ----
dir.create("Figure", showWarnings = FALSE)
ggsave("Figure/IL1R2_spline_vertical.svg", final_plot, width = 10, height = 6, dpi = 300)
cat("✅ Figure saved: Figure/IL1R2_spline_vertical.svg\n")
