## ===============================
## MS1 tertiles → Table 1
## ===============================

rm(list = ls())
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2); library(tableone)
})

# 1) Read MS1 % from CIBERSORTx
ms1_df <- readr::read_csv("Original_data/CIBERSORTx_percent_MS1.csv",
                          show_col_types = FALSE)
# In your earlier code it's percent (0–100). Make sure it's numeric:
ms1_df <- ms1_df %>%
  rename(Sample = Mixture) %>%
  mutate(MS1 = as.numeric(MS1))

stopifnot(!anyNA(ms1_df$MS1))

# 2) Map Sample → ID
dfwb_sel <- readr::read_csv("Original_data/dfwb_sel.csv", show_col_types = FALSE)
ms1_map <- ms1_df %>%
  dplyr::left_join(
    dfwb_sel %>% dplyr::select(Sample, ID),
    by = "Sample"
  )


# Optional sanity check
if (any(is.na(ms1_map$ID))) {
  message("Warning: some Samples in CIBERSORTx output are missing in dfwb_sel. They will be dropped.")
}
ms1_map <- ms1_map %>% filter(!is.na(ID))

# 3) Create MS1 tertiles by rank (equal sized groups)
#    → 1 = Low, 2 = Mid, 3 = High
ms1_map <- ms1_map %>%
  mutate(
    MS1_group_num = dplyr::ntile(MS1, 3),
    MS1_group = factor(MS1_group_num,
                       levels = c(1,2,3),
                       labels = c("Low", "Mid", "High"))
  )

# Quick counts
print(table(ms1_map$MS1_group))

# 4) Merge into studydata
studydata <- readr::read_csv("Original_data/studydata.csv", show_col_types = FALSE)

studydata_with_ms1 <- studydata %>%
  mutate(EB_id = as.character(EB_id)) %>%
  inner_join(ms1_map %>% transmute(EB_id = as.character(ID),
                                   Sample,
                                   MS1, MS1_group),
             by = "EB_id")

# Export the merged clinical+MS1 file
write.csv(studydata_with_ms1,
          file.path("Original_data", "studydata_with_ms1.csv"),
          row.names = FALSE)
cat("Saved: Original_data/studydata_with_ms1.csv\n")


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
             "hospdeath",	"mortality_d30", "mortality_d90","MS1_group")

catvars <- c("group","flow_Data", "spectral", "inclusion_hospital", "gender",	"ethnic_group#White",	
             "hypertension",	"cpd",	"COPD",	"diabetes",	
             "ccd","ckd",	"mneoplasm",	"immune_sup",	"cnd",	
             "Oxygen_therapy_1", 
             "antibiotic_seven_days", 
             "ICU_Medium_Care_admission_1",	 
             "Non_invasive_ventilation_1",	"Invasive_ventilation_1",	
             "bacteremia",
             "hospdeath",	"mortality_d30", "mortality_d90", "pathogen_cultured","MS1_group")

nonnormal <- c("sampling_time",		"age_yrs", "BMI", "symptom_days",  "CCI", "qSOFA_score",	"MEWS_score",	"CURB_score",	"PSI_new",
               "crp_1_1",	 "WBC_2_1",	"Neutrophil_unit_1",	"Lymphocyte_1_1",	
               "Creatinine_value_1",
               "Platelets_value_1",	
               "Blood_Urea_Nitrogen_value_1",	"length_of_oxygen","length_of_stay" ,"icu_stay", "lenght_of_intubation")

tab1     <- CreateTableOne(
  vars        = allvars, 
  data        = studydata_with_ms1, 
  strata = "MS1_group",
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
export_file_name <- "table1_ms1.csv" 
write.csv(tab1_df, file.path(destination_folder, export_file_name), row.names = FALSE)
cat("File exported to:", file.path(destination_folder, export_file_name), "\n")


## =========================================================
## MS1 version – keep calculation, only change facet titles
## =========================================================

rm(list = ls())

## ===== Packages =====
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(stringr); library(forcats)
  library(ggplot2); library(rstatix); library(purrr); library(patchwork)
})

## ===== Data input =====
studydata_with_ms1 <- read_csv("Original_data/studydata_with_ms1.csv", show_col_types = FALSE)
whole_blood_stimuli_long <- read_csv("Original_data/whole_blood_stimuli_long.csv", show_col_types = FALSE)
if ("HiIL6" %in% names(whole_blood_stimuli_long)) whole_blood_stimuli_long$HiIL6 <- NULL

## ===== ID key =====
id_key <- studydata_with_ms1 %>%
  transmute(
    EB_id     = as.character(EB_id),
    MS1       = as.numeric(MS1),
    MS1_group = factor(MS1_group, levels = c("Low","Mid","High"))
  )

# 如果没有分组则自动三分位
if (any(is.na(id_key$MS1_group))) {
  id_key <- id_key %>%
    mutate(
      MS1_group_num = ntile(MS1, 3),
      MS1_group = factor(MS1_group_num, levels = c(1,2,3), labels = c("Low","Mid","High"))
    ) %>% select(-MS1_group_num)
}

## ===== Join & clean =====
wb <- whole_blood_stimuli_long %>%
  mutate(ID = as.character(ID)) %>%
  semi_join(id_key, by = c("ID" = "EB_id")) %>%
  left_join(id_key,  by = c("ID" = "EB_id"))

stim_keep <- c("KP","PP","LPS")
group_var <- "MS1_group"
analysis_scale <- "log10"

wb_clean <- wb %>%
  filter(!is.na(marker), !is.na(stimulation), !is.na(.data[[group_var]])) %>%
  filter(stimulation %in% stim_keep) %>%
  mutate(
    value_raw = as.character(value),
    v_dot   = parse_number(value_raw, locale = locale(decimal_mark = ".", grouping_mark = ",")),
    v_comma = parse_number(value_raw, locale = locale(decimal_mark = ",", grouping_mark = ".")),
    value   = coalesce(v_dot, v_comma),
    value   = if_else(str_detect(value_raw, "^\\s*<\\s*\\d+(?:[\\.,]\\d+)?\\s*$"),
                      parse_number(str_replace(value_raw, "<","")), value),
    value   = if_else(str_detect(value_raw, "^\\s*>\\s*\\d+(?:[\\.,]\\d+)?\\s*$"),
                      parse_number(str_replace(value_raw, ">","")), value),
    value_num = if (analysis_scale == "log10") log10(pmax(value, 1e-9)) else value,
    stimulation = factor(stimulation, levels = stim_keep),
    !!group_var := factor(.data[[group_var]], levels = c("Low","Mid","High"))
  ) %>%
  select(stimulation, marker, value_num, !!sym(group_var))

## ===== 计算部分（保持不变） =====
pairs_df <- tribble(
  ~first, ~second,
  "Low",  "Mid",
  "Mid",  "High",
  "Low",  "High"
)

effsize_df <- purrr::pmap_dfr(pairs_df, function(first, second) {
  dat <- wb_clean %>%
    filter(.data[[group_var]] %in% c(first, second)) %>%
    mutate(!!group_var := fct_relevel(.data[[group_var]], second, first))
  
  es <- dat %>%
    group_by(stimulation, marker) %>%
    rstatix::cohens_d(reformulate(group_var, response = "value_num"),
                      var.equal = FALSE, hedges.correction = TRUE, ci = TRUE) %>%
    ungroup() %>% rename(g = effsize)
  
  tt <- dat %>%
    group_by(stimulation, marker) %>%
    rstatix::t_test(reformulate(group_var, response = "value_num"), var.equal = FALSE) %>%
    ungroup() %>%
    mutate(p_BH = p.adjust(p, method = "BH")) %>%
    select(stimulation, marker, statistic, p_BH)
  
  es %>%
    left_join(tt, by = c("stimulation","marker")) %>%
    mutate(
      contrast_lab = sprintf("MS1-%s vs MS1-%s\n(positive = higher in first group)", first, second),
      p_BH.signif = case_when(
        is.na(p_BH) ~ "ns",
        p_BH < 1e-4 ~ "****",
        p_BH < 1e-3 ~ "***",
        p_BH < 1e-2 ~ "**",
        p_BH < 5e-2 ~ "*",
        TRUE ~ "ns"
      ),
      dir = case_when(statistic > 0 ~ "UP",
                      statistic < 0 ~ "DOWN",
                      TRUE ~ NA_character_),
      signif_dir = ifelse(p_BH.signif == "ns" | is.na(dir), "ns", paste(p_BH.signif, dir))
    )
})

## ===== 改 facet 标注（只改标签，不改计算） =====
lab_map <- c(
  "MS1-Low vs MS1-Mid\n(positive = higher in first group)" = "MS1-Mid vs MS1-Low\n(positive = higher in first group)",
  "MS1-Mid vs MS1-High\n(positive = higher in first group)" = "MS1-High vs MS1-Mid\n(positive = higher in first group)",
  "MS1-Low vs MS1-High\n(positive = higher in first group)" = "MS1-High vs MS1-Low\n(positive = higher in first group)"
)

stim_labels <- c(KP = "Kleb.", PP = "Pneumo.", LPS = "LPS")
marker_pretty <- c(
  "IL_1beta"="IL-1β","IL_1RA"="IL-1RA","IL_8_1"="IL-8","IL_6"="IL-6",
  "IL_10"="IL-10","TNF"="TNF","CCL2"="CCL2","CCL3"="CCL3","CCL4"="CCL4"
)

marker_order <- c("CCL2","CCL3","CCL4","IL-1β","IL-6","TNF","IL-8","IL-1RA","IL-10")
contrast_order <- c(
  "MS1-Low vs MS1-Mid\n(positive = higher in first group)",
  "MS1-Mid vs MS1-High\n(positive = higher in first group)",
  "MS1-Low vs MS1-High\n(positive = higher in first group)"
)

effsize_df <- effsize_df %>%
  mutate(
    stimulation = fct_relabel(stimulation, ~ stim_labels[.x] %||% .x),
    marker = factor(marker_pretty[as.character(marker)] %||% as.character(marker),
                    levels = marker_order),
    contrast_lab = factor(contrast_lab, levels = contrast_order)
  )

## ===== 绘图 =====
color_mapping <- c(
  "ns"="#888888",
  "* UP"="#ffd4c4","** UP"="#FF9F94","*** UP"="#ff5c33","**** UP"="#cc2900",
  "* DOWN"="#c2d1ff","** DOWN"="#99b3ff","*** DOWN"="#668cff","**** DOWN"="blue"
)

strip_bg <- effsize_df %>%
  group_by(stimulation, contrast_lab) %>%
  distinct(marker) %>%
  arrange(stimulation, contrast_lab, marker) %>%
  mutate(row_id = row_number(),
         ymin = row_id - 0.5, ymax = row_id + 0.5) %>%
  filter(row_id %% 2 == 1) %>%
  ungroup()

p_all <- ggplot(effsize_df, aes(x = g, y = marker)) +
  geom_rect(data = strip_bg, inherit.aes = FALSE,
            aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
            fill = "gray95", colour = NA) +
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
  facet_grid(rows = vars(stimulation),
             cols = vars(contrast_lab),
             scales = "free_x",
             labeller = labeller(contrast_lab = as_labeller(lab_map))) +
  labs(
    x = "Hedges’ g (positive = higher in first group), 95% CI",
    y = NULL,
    caption = "Positive values indicate higher cytokine response in the first group (the group before 'vs')."
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y   = element_text(face = "bold"),
    strip.text.x  = element_text(face = "bold"),
    strip.text.y  = element_text(face = "bold"),
    panel.grid    = element_blank(),
    axis.line.x   = element_line(size = 0.25),
    legend.position = "none",
    plot.margin   = margin(10, 30, 10, 10),
    panel.spacing.x = unit(10, "pt"),
    panel.spacing.y = unit(8, "pt"),
    plot.caption = element_text(hjust = 0, size = 10, color = "gray30")
  )

print(p_all)
ggsave("Figure/ms1_hedges_labelReversed.svg", p_all, width = 12, height = 8, dpi = 300)

## ===== Box plot (保持原逻辑) =====
stim_now <- "KP"
dat_kleb <- wb_clean %>%
  filter(stimulation == stim_now) %>%
  mutate(marker = factor(marker),
         MS1_group = as.factor(MS1_group))

colors <- c("Low"="#f28e8e", "Mid"="#5bc0de", "High"="#a1d99b")

make_plot <- function(first, second, contrast_lab) {
  dat_sub <- dat_kleb %>%
    filter(MS1_group %in% c(first, second)) %>%
    mutate(MS1_group = fct_relevel(MS1_group, first, second))
  
  ggplot(dat_sub, aes(x = MS1_group, y = value_num, fill = MS1_group)) +
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

p1 <- make_plot("Low", "Mid",  "MS1-Mid vs MS1-Low")
p2 <- make_plot("Mid", "High", "MS1-High vs MS1-Mid")
p3 <- make_plot("Low", "High", "MS1-High vs MS1-Low")

final_plot <- p1 / p2 / p3 +
  plot_annotation(
    title = "Kleb. stimulation — MS1 contrasts across all markers",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  )

print(final_plot)
ggsave("Figure/Box_allmarkers_MS1_labelReversed.svg", final_plot, width = 12, height = 8, dpi = 300)




#for the continous score####
## ===============================
##  MS1 spline plots (lm + rcs)
## ===============================

library(dplyr)
library(ggplot2)
library(rms)
library(readr)
library(patchwork)

# ---- Load data ----
studydata_with_ms1 <- read_csv("Original_data/studydata_with_ms1.csv", show_col_types = FALSE)
whole_blood_stimuli_long <- read_csv("Original_data/whole_blood_stimuli_long.csv", show_col_types = FALSE)

# ---- Markers and stimulations of interest ----
markers_keep <- c("TNF", "IL_1RA")
stim_keep <- c("LPS", "PP", "KP")

# ---- Prepare dataset ----
wb <- whole_blood_stimuli_long %>%
  filter(marker %in% markers_keep, stimulation %in% stim_keep) %>%
  mutate(ID = as.character(ID)) %>%
  left_join(
    studydata_with_ms1 %>% 
      transmute(EB_id = as.character(EB_id), MS1 = as.numeric(MS1)),
    by = c("ID" = "EB_id")
  ) %>%
  # convert to character before parsing (fixes your error)
  mutate(
    value_raw = as.character(value),
    value     = readr::parse_number(value_raw),
    log10_value = log10(pmax(value, 1e-9))
  ) %>%
  filter(!is.na(MS1), !is.na(log10_value))

# ---- Simple spline plot function ----
plot_spline <- function(stim_now, marker_now) {
  dat <- wb %>% filter(stimulation == stim_now, marker == marker_now)
  
  # fit spline model (lm + rcs)
  fit <- lm(log10_value ~ rcs(MS1, 3), data = dat)
  
  # predict smooth line
  newdat <- data.frame(MS1 = seq(min(dat$MS1), max(dat$MS1), length.out = 200))
  p <- predict(fit, newdata = newdat, se.fit = TRUE)
  newdat$fit <- p$fit
  newdat$upper <- p$fit + 1.96 * p$se.fit
  newdat$lower <- p$fit - 1.96 * p$se.fit
  
  # draw
  ggplot(newdat, aes(x = MS1, y = fit)) +
    geom_line(color = "#0072B2", linewidth = 1) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#0072B2", alpha = 0.2) +
    labs(
      title = paste(stim_now, "—", marker_now),
      x = "MS1 score",
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

# ---- Combine in one figure ----
# Column 1 = TNF (LPS / PP / KP top→bottom)
# Column 2 = IL-1RA (LPS / PP / KP top→bottom)
final_plot <- (plots[[1]] / plots[[2]] / plots[[3]]) | 
  (plots[[4]] / plots[[5]] / plots[[6]])


print(final_plot)
ggsave("Figure/MS1_spline.svg", final_plot, width = 10, height = 6, dpi = 300)

