rm(list = ls ())

## =========================
## 0) Load dependencies (must be installed already)
## =========================
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
               tableone, kmed, survival, survminer, skimr, lattice, reshape2, grDevices, ggrepel)


SRSq <- import("Original_data/SRSq_and_SRS12_only.csv") 

SRSq7 <- SRSq %>%
  select(ID, matches("^SRS.*_7$"))

studydata <- readr::read_csv("Original_data/studydata.csv", show_col_types = FALSE)
whole_blood_stimuli_long <- readr::read_csv("Original_data/whole_blood_stimuli_long.csv", show_col_types = FALSE)
if ("HiIL6" %in% names(whole_blood_stimuli_long)) whole_blood_stimuli_long$HiIL6 <- NULL

## ---------- 1) Clean + make SRSq_7 tertiles ----------
srs_key <- SRSq7 %>%
  transmute(
    EB_id      = as.character(ID),
    SRS_label  = factor(SRS_7, levels = c("SRS1","SRS2")),
    SRSq_7     = as.numeric(SRSq_7)
  ) %>%
  mutate(
    SRSq7_group_num = ntile(SRSq_7, 3),
    SRSq7_group     = factor(SRSq7_group_num, levels = c(1,2,3),
                             labels = c("Low","Mid","High"))
  ) %>% select(-SRSq7_group_num)

## ---------- 2) Merge to clinical data ----------
studydata_srs <- studydata %>%
  mutate(EB_id = as.character(EB_id)) %>%
  inner_join(srs_key, by = "EB_id")

write.csv(studydata_srs, "Original_data/studydata_with_SRSq7.csv", row.names = FALSE)
cat("Saved: Original_data/studydata_with_SRSq7.csv\n")

## ---------- 3) Table 1 by SRSq7 tertiles ----------
allvars <- c("group","flow_Data", "spectral", "inclusion_hospital", "sampling_time",
             "age_yrs","gender","ethnic_group#White","BMI","symptom_days",
             "CCI","hypertension","cpd","COPD","diabetes",
             "ccd","ckd","mneoplasm","immune_sup","cnd",
             "qSOFA_score","MEWS_score","CURB_score","PSI_new",
             "crp_1_1","WBC_2_1","Neutrophil_unit_1","Lymphocyte_1_1",
             "Creatinine_value_1","Platelets_value_1","Blood_Urea_Nitrogen_value_1",
             "Oxygen_therapy_1","length_of_oxygen","antibiotic_seven_days",
             "length_of_stay","ICU_Medium_Care_admission_1","icu_stay",
             "Non_invasive_ventilation_1","Invasive_ventilation_1","lenght_of_intubation",
             "bacteremia","pathogen_cultured","hospdeath","mortality_d30","mortality_d90",
             "SRSq7_group")

catvars <- c("group","flow_Data","spectral","inclusion_hospital","gender","ethnic_group#White",
             "hypertension","cpd","COPD","diabetes","ccd","ckd","mneoplasm","immune_sup","cnd",
             "Oxygen_therapy_1","antibiotic_seven_days","ICU_Medium_Care_admission_1",
             "Non_invasive_ventilation_1","Invasive_ventilation_1","bacteremia",
             "hospdeath","mortality_d30","mortality_d90","pathogen_cultured","SRSq7_group")

nonnormal <- c("sampling_time","age_yrs","BMI","symptom_days","CCI","qSOFA_score","MEWS_score",
               "CURB_score","PSI_new","crp_1_1","WBC_2_1","Neutrophil_unit_1","Lymphocyte_1_1",
               "Creatinine_value_1","Platelets_value_1","Blood_Urea_Nitrogen_value_1",
               "length_of_oxygen","length_of_stay","icu_stay","lenght_of_intubation")

tab1 <- CreateTableOne(
  vars = allvars, data = studydata_srs,
  strata = "SRSq7_group", factorVars = catvars, test = TRUE
)

print(tab1, nonnormal = nonnormal, quote = FALSE, noSpaces = TRUE, smd = TRUE, missing = TRUE)
tab1_df <- as.data.frame(print(tab1, nonnormal = nonnormal, quote = FALSE, noSpaces = TRUE, smd = TRUE, missing = TRUE))
tab1_df$tablename <- rownames(tab1_df)
tab1_df <- tab1_df[, c("tablename", setdiff(names(tab1_df), "tablename"))]
write.csv(tab1_df, "Original_data/table1_SRSq7.csv", row.names = FALSE)
cat("Saved: Original_data/table1_SRSq7.csv\n")

## ---------- 4) Join markers to SRSq7 tertiles ----------
id_key <- srs_key %>% select(EB_id, SRSq_7, SRSq7_group, SRS_label)

wb <- whole_blood_stimuli_long %>%
  mutate(ID = as.character(ID)) %>%
  semi_join(id_key, by = c("ID" = "EB_id")) %>%
  left_join(id_key,  by = c("ID" = "EB_id"))

## ---------- 5) Effect sizes across tertiles (Hedges’ g) ----------
stim_keep <- c("KP","PP","LPS")           # keep your 3 stimuli
group_var <- "SRSq7_group"
analysis_scale <- "log10"

wb_clean <- wb %>%
  filter(stimulation %in% stim_keep, !is.na(marker), !is.na(.data[[group_var]])) %>%
  mutate(
    value_num = log10(pmax(readr::parse_number(as.character(value)), 1e-9)),
    stimulation = factor(stimulation, levels = stim_keep),
    !!group_var := factor(.data[[group_var]], levels = c("Low","Mid","High"))
  ) %>%
  select(stimulation, marker, value_num, !!sym(group_var))

# pairwise contrasts
lvls <- levels(wb_clean[[group_var]])
contrast_df <- as.data.frame(t(combn(lvls, 2))); names(contrast_df) <- c("ref","comp")

# (optional) ensure plyr is not masking dplyr
if ("package:plyr" %in% search()) detach("package:plyr", unload = TRUE)

effsize_df <- purrr::pmap_dfr(contrast_df, function(ref, comp){
  dat <- wb_clean %>%
    dplyr::filter(.data[[group_var]] %in% c(ref, comp)) %>%
    dplyr::mutate(!!group_var := forcats::fct_relevel(.data[[group_var]], ref, comp))
  
  es <- dat %>%
    dplyr::group_by(stimulation, marker) %>%
    rstatix::cohens_d(reformulate(group_var, response = "value_num"),
                      var.equal = FALSE, hedges.correction = TRUE, ci = TRUE) %>%
    dplyr::ungroup() %>%
    dplyr::rename(g = effsize)     # <- use dplyr::rename
  
  tt <- dat %>%
    dplyr::group_by(stimulation, marker) %>%
    rstatix::t_test(reformulate(group_var, response = "value_num"), var.equal = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(p_BH = p.adjust(p, method = "BH")) %>%
    dplyr::select(stimulation, marker, statistic, p_BH)
  
  es %>% dplyr::left_join(tt, by = c("stimulation","marker")) %>%
    dplyr::mutate(
      contrast_lab = sprintf("%s vs %s (ref = %s)", comp, ref, ref),
      p_BH.signif = dplyr::case_when(
        p_BH < 1e-4 ~ "****", p_BH < 1e-3 ~ "***", p_BH < 1e-2 ~ "**",
        p_BH < 5e-2 ~ "*", TRUE ~ "ns"
      ),
      dir = dplyr::case_when(statistic > 0 ~ "UP", statistic < 0 ~ "DOWN", TRUE ~ NA_character_),
      signif_dir = ifelse(p_BH.signif == "ns" | is.na(dir), "ns", paste(p_BH.signif, dir))
    )
})


# order panels: High vs Low first
contrast_order <- c(
  "High vs Low (ref = Low)", "Mid vs Low (ref = Low)", "High vs Mid (ref = Mid)"
)
effsize_df <- effsize_df %>%
  mutate(
    contrast_lab = factor(contrast_lab,
                          levels = contrast_order)
  )

# striped background
strip_bg <- effsize_df %>%
  group_by(stimulation, contrast_lab) %>%
  distinct(marker) %>%
  arrange(stimulation, contrast_lab, marker) %>%
  mutate(row_id = row_number(), ymin = row_id - 0.5, ymax = row_id + 0.5) %>%
  filter(row_id %% 2 == 1) %>% ungroup()

color_mapping <- c(
  "ns"="#888888",
  "* UP"="#ffd4c4","** UP"="#FF9F94","*** UP"="#ff5c33","**** UP"="#cc2900",
  "* DOWN"="#c2d1ff","** DOWN"="#99b3ff","*** DOWN"="#668cff","**** DOWN"="blue"
)

p_srsq_hedges <- ggplot(effsize_df, aes(x = g, y = marker)) +
  geom_rect(data = strip_bg, inherit.aes = FALSE,
            aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
            fill = "gray95", colour = NA) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_segment(aes(x = conf.low, xend = conf.high, yend = marker), linewidth = 0.4) +
  geom_point(aes(fill = signif_dir), shape = 22, size = 3) +
  scale_fill_manual(values = color_mapping, na.value = "#888888") +
  facet_grid(rows = vars(stimulation), cols = vars(contrast_lab), scales = "free_x") +
  labs(x = "Hedges’ g, 95% CI", y = NULL,
       title = "Stimulated whole blood: SRSq_7 tertiles (effect sizes)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        panel.grid = element_blank(),
        legend.position = "none")

dir.create("Figure", showWarnings = FALSE)
ggsave("Figure/SRSq7_hedges.svg", p_srsq_hedges, width = 12, height = 8, dpi = 300)

## ---------- 6) Primary boxplots by tertiles ----------
wb2 <- wb %>%
  mutate(
    stimulation = factor(stimulation,
                         levels = c("M","KP","PP","LPS"),
                         labels = c("Medium","Kleb.","Pneumo.","LPS")),
    SRSq7_group = factor(SRSq7_group, levels = c("Low","Mid","High"))
  )

p_srsq_box <- ggplot(wb2, aes(x = stimulation, y = value, fill = SRSq7_group)) +
  geom_boxplot(position = position_dodge(width = 0.75), width = 0.6,
               alpha = 0.55, colour = "#4a4a4a", outlier.shape = NA) +
  geom_point(aes(color = SRSq7_group),
             position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
             size = 1.4, alpha = 0.65) +
  scale_y_continuous(trans = "log10") +
  facet_wrap(~ marker, scales = "free_y", ncol = 3) +
  labs(x = NULL, y = "Cytokine (log10 scale)",
       title = "Whole-blood stimulation: SRSq_7 tertiles") +
  theme_classic(base_size = 12) +
  theme(strip.text = element_text(face = "bold"),
        legend.position = "top") +
  guides(color = "none")

ggsave("Figure/SRSq7_boxplots.svg", p_srsq_box, width = 12, height = 8, dpi = 300)

## ---------- 7) Significance grid (Wilcoxon, BH) ----------
wb3 <- wb2
my_comps <- list(c("Low","Mid"), c("Low","High"), c("Mid","High"))

stat_df <- wb3 %>%
  group_by(marker, stimulation) %>%
  group_modify(~{
    d <- .x
    y_max <- max(d$value, na.rm = TRUE)
    rstatix::wilcox_test(d, value ~ SRSq7_group, comparisons = my_comps, p.adjust.method = "BH") %>%
      rstatix::add_significance(p.col = "p.adj") %>%
      mutate(
        pair = paste(group1, group2, sep = "-"),
        rank = match(pair, c("Low-Mid","Low-High","Mid-High")),
        y.position = y_max * (1.15 + 0.08 * (rank - 1))
      )
  }) %>% ungroup()

p_srsq_sig <- ggplot(wb3, aes(x = SRSq7_group, y = value, fill = SRSq7_group)) +
  geom_boxplot(width = 0.65, alpha = 0.55, colour = "#4a4a4a", outlier.shape = NA) +
  geom_jitter(aes(color = SRSq7_group), width = 0.12, size = 1.2, alpha = 0.6) +
  scale_y_continuous(trans = "log10") +
  facet_grid(marker ~ stimulation, scales = "free_y") +
  labs(x = NULL, y = "Cytokine (log10 scale)",
       title = "Whole-blood stimulation by SRSq_7 tertiles (BH-adjusted)") +
  theme_classic(base_size = 12) +
  theme(strip.text = element_text(face = "bold"),
        legend.position = "top") +
  guides(color = "none") +
  ggpubr::stat_pvalue_manual(
    stat_df, label = "p.adj.signif",
    hide.ns = TRUE, tip.length = 0.01, bracket.size = 0.3, size = 3
  )

ggsave("Figure/SRSq7_significance_grid.svg", p_srsq_sig, width = 12, height = 10, dpi = 300)

cat("\nSaved:\n- Original_data/table1_SRSq7.csv\n- Figure/SRSq7_hedges.svg\n- Figure/SRSq7_boxplots.svg\n- Figure/SRSq7_significance_grid.svg\n")
