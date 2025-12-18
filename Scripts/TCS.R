# ==== Time to Clinical Stability (TTCS) ====
rm(list = ls())

## install.packages("pacman")
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
               tableone, kmed, survival, survminer, skimr, lattice, reshape2, grDevices, ggrepel)

install.packages("dplyr")
library(dplyr)
library(MASS)
library(magrittr)
library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)
library(readxl)

# ---- Load packages ----
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, rio, lubridate, janitor, openxlsx, dplyr, stringr, readxl)

# ---- 1) Read study data ----
# studydata should contain: EB_id, admission_dt1, discharge_date_local, day28_pat_status, date_of_death
studydata <- import("Original_data/studydata.csv") %>%
  mutate(
    admission_dt1        = suppressWarnings(ymd(admission_dt1)),
    discharge_date_local = suppressWarnings(ymd(discharge_date_local)),
    date_of_death        = suppressWarnings(ymd(date_of_death))
  ) %>%
  filter(group == "CAP")
ids_keep <- unique(studydata$EB_id)

# ---- 2) Read Excel sheet ----
xlsx_path  <- "Original_data/ELDER-BIOME_excel_export_20230418012733.xlsx"
all_sheets <- readxl::excel_sheets(xlsx_path)
sheet_name <- "Clinical_stabilty_check"
if (!sheet_name %in% all_sheets) {
  stop("Sheet not found: ", sheet_name, "\nAvailable sheets:\n", paste(all_sheets, collapse = ", "))
}

df_tcs <- readxl::read_excel(xlsx_path, sheet = sheet_name, col_types = "text") %>%
  clean_names() %>%
  dplyr::mutate(clin_stab_hosp = as.numeric(clin_stab_hosp)) %>%
  dplyr::filter(participant_id %in% ids_keep) %>%
  dplyr::mutate(record_id = as.character(participant_id) |> trimws())

# ---- 3) Merge with studydata ----
studydata_keyed <- studydata %>%
  dplyr::mutate(EB_id = as.character(EB_id) |> trimws()) %>%
  dplyr::select(EB_id, admission_dt1, discharge_date_local, day28_pat_status, date_of_death) %>%
  dplyr::rename(record_id = EB_id)

df_tcs <- left_join(df_tcs, studydata_keyed, by = "record_id") %>%
  mutate(
    clin_stab_date = suppressWarnings(dmy(clin_stab_date)),
    los = as.numeric(discharge_date_local - admission_dt1),
    los = ifelse(is.na(los), NA_real_, pmax(los, 0)),
    day28_pat_status_num = suppressWarnings(as.numeric(day28_pat_status)),
    day28_pat_status_num = ifelse(
      is.na(day28_pat_status_num) &
        str_detect(tolower(as.character(day28_pat_status)), "deceased|dead|死亡|died"),
      5, day28_pat_status_num
    )
  )

# ---- 4) Define stability flags (Halm criteria) ----
make_stability_flags <- function(df, temp_cutoff = 37.2) {
  df %>%
    mutate(
      clin_stab_temp   = suppressWarnings(as.numeric(clin_stab_temp)),
      clin_stab_hr     = suppressWarnings(as.numeric(clin_stab_hr)),
      clin_stab_resp   = suppressWarnings(as.numeric(clin_stab_resp)),
      clin_stab_sbp    = suppressWarnings(as.numeric(clin_stab_sbp)),
      clin_stab_oxygen = suppressWarnings(as.numeric(clin_stab_oxygen)),
      clin_stab_add_oxygen_std = tolower(trimws(as.character(clin_stab_add_oxygen))),
      cs_temp = ifelse(is.na(clin_stab_temp) | clin_stab_temp <= temp_cutoff, 1L, 0L),
      cs_hr   = ifelse(is.na(clin_stab_hr)   | clin_stab_hr   <= 100, 1L, 0L),
      cs_rr   = ifelse(is.na(clin_stab_resp) | clin_stab_resp <= 24, 1L, 0L),
      cs_bp   = ifelse(is.na(clin_stab_sbp)  | clin_stab_sbp  >= 90, 1L, 0L),
      cs_ox   = ifelse(clin_stab_add_oxygen_std == "no" & clin_stab_oxygen >= 90, 1L, 0L),
      stable  = ifelse(cs_temp==1L & cs_hr==1L & cs_rr==1L & cs_bp==1L & cs_ox==1L, 1L, 0L),
      stable  = ifelse(!is.na(clin_stab_hosp) & clin_stab_hosp == 0, 0L, stable)
    )
}

df_tcs_old <- make_stability_flags(df_tcs, temp_cutoff = 37.2)
df_tcs_new <- make_stability_flags(df_tcs, temp_cutoff = 37.8)

# ---- 5) Pick best record per patient ----
pick_one_row <- function(df, out_name = "ttcs") {
  df %>%
    arrange(record_id, desc(stable), is.na(clin_stab_hosp), clin_stab_hosp) %>%
    distinct(record_id, .keep_all = TRUE) %>%
    transmute(
      record_id,
      !!out_name := case_when(
        clin_stab_hosp == 0 & stable == 0 & day28_pat_status_num == 5 ~ 29,
        clin_stab_hosp == 0 & stable == 0 ~ pmin(as.numeric(los), 28),
        TRUE ~ pmin(as.numeric(clin_stab_hosp), 28)
      )
    )
}

df_stable_old <- pick_one_row(df_tcs_old, "ttcs_halm_372")
df_stable_new <- pick_one_row(df_tcs_new, "ttcs_halm_378")

# ---- 6) Final TTCS table ----
studydata_key <- studydata %>%
  mutate(EB_id = as.character(EB_id) |> trimws()) %>%
  dplyr::select(EB_id, admission_dt1, discharge_date_local, day28_pat_status)

ttcs <- full_join(df_stable_old, df_stable_new, by = "record_id") %>%
  dplyr::rename(EB_id = record_id) %>%
  dplyr::left_join(studydata_key, by = "EB_id") %>%
  dplyr::mutate(
    ttcs_halm_372_days  = ttcs_halm_372,
    ttcs_halm_378_days  = ttcs_halm_378,
    ttcs_halm_372_hours = round(ttcs_halm_372 * 24, 2),
    ttcs_halm_378_hours = round(ttcs_halm_378 * 24, 2)
  ) %>%
  dplyr::relocate(EB_id,
           ttcs_halm_372, ttcs_halm_378,
           ttcs_halm_372_days, ttcs_halm_378_days,
           ttcs_halm_372_hours, ttcs_halm_378_hours)


ttcs <- ttcs %>%
  mutate(
    ttcs_group = case_when(
      !is.na(ttcs_halm_372_days) & ttcs_halm_372_days <= 3 ~ "short",
      !is.na(ttcs_halm_372_days) & ttcs_halm_372_days > 3  ~ "long",
      TRUE ~ NA_character_
    )
  )

ttcs <- ttcs %>%
  mutate(ttcs_group = recode(ttcs_group,
                             "short" = "≤3 days",
                             "long"  = ">3 days"))

ggplot(ttcs, aes(x = ttcs_halm_372_days, fill = ttcs_group)) +
  geom_histogram(binwidth = 2, color = "black", alpha = 0.7, position = "identity") +
  scale_fill_manual(values = c("≤3 days" = "#1f77b4", ">3 days" = "#ff7f0e")) +
  labs(
    title = "Time to Clinical Stability",
    x = "TTCS days",
    y = "Frequency",
    fill = "Group"
  ) +
  theme_classic(base_size = 14)



#merge
merged_data <- ttcs %>%
  filter(!is.na(ttcs_halm_372_days)) %>%
  mutate(EB_id = as.character(EB_id)) %>%
  left_join(studydata %>% mutate(EB_id = as.character(EB_id)), by = "EB_id")


## table1
## analyse errors 
allvars <- c("ttcs_group","flow_Data", "spectral", "inclusion_hospital", "sampling_time",		"age_yrs",	"gender",	"ethnic_group#White",	
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
             "hospdeath",	"mortality_d30", "mortality_d90")

catvars <- c("ttcs_group","flow_Data", "spectral", "inclusion_hospital", "gender",	"ethnic_group#White",	
             "hypertension",	"cpd",	"COPD",	"diabetes",	
             "ccd","ckd",	"mneoplasm",	"immune_sup",	"cnd",	
             "Oxygen_therapy_1", 
             "antibiotic_seven_days", 
             "ICU_Medium_Care_admission_1",	 
             "Non_invasive_ventilation_1",	"Invasive_ventilation_1",	
             "bacteremia",
             "hospdeath",	"mortality_d30", "mortality_d90", "pathogen_cultured")

nonnormal <- c("sampling_time",		"age_yrs", "BMI", "symptom_days",  "CCI", "qSOFA_score",	"MEWS_score",	"CURB_score",	"PSI_new",
               "crp_1_1",	 "WBC_2_1",	"Neutrophil_unit_1",	"Lymphocyte_1_1",	
               "Creatinine_value_1",
               "Platelets_value_1",	
               "Blood_Urea_Nitrogen_value_1",	"length_of_oxygen","length_of_stay" ,"icu_stay", "lenght_of_intubation")

tab1     <- CreateTableOne(
  vars        = allvars, 
  data        = merged_data, 
  strata = "ttcs_group",
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
export_file_name <- "table1_tcs.csv" 
write.csv(tab1_df, file.path(destination_folder, export_file_name), row.names = FALSE)
cat("File exported to:", file.path(destination_folder, export_file_name), "\n")



#hedge and vilon plot
## ===== 0) Packages =====
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, rstatix, ggplot2, ggrepel)

## If plyr is loaded, suggest unloading it to avoid conflicts with dplyr
if ("package:plyr" %in% search()) detach("package:plyr", unload = TRUE)

## ===== 1) Parameters =====
stim_keep <- c("KP","PP","LPS")

## ===== 2) Clean luminex_long$value into numeric =====
# Notes: Compatible with decimal comma, text with symbols (e.g. "<5", ">100", "12,3"); thresholds are treated as numeric
luminex_long <- luminex_long %>%
  mutate(
    value_raw = as.character(value),
    value_num_dot = readr::parse_number(
      value_raw,
      locale = readr::locale(decimal_mark = ".", grouping_mark = ",")
    ),
    value_num_comma = readr::parse_number(
      value_raw,
      locale = readr::locale(decimal_mark = ",", grouping_mark = ".")
    ),
    value = dplyr::coalesce(value_num_dot, value_num_comma),
    value = dplyr::if_else(
      stringr::str_detect(value_raw, "^\\s*<\\s*\\d+(?:[\\.,]\\d+)?\\s*$"),
      readr::parse_number(stringr::str_replace(value_raw, "<", "")),
      value
    ),
    value = dplyr::if_else(
      stringr::str_detect(value_raw, "^\\s*>\\s*\\d+(?:[\\.,]\\d+)?\\s*$"),
      readr::parse_number(stringr::str_replace(value_raw, ">", "")),
      value
    )
  ) %>%
  select(-value_num_dot, -value_num_comma)

## (Optional) Check raw entries that are still NA
# luminex_long %>% filter(is.na(value) & !is.na(value_raw)) %>% count(marker)

## ===== 3) Merge merged_data, and only keep KP/PP/LPS =====
dat <- luminex_long %>%
  filter(stimulation %in% stim_keep) %>%
  inner_join(
    merged_data %>%
      transmute(ID = as.character(EB_id),
                ttcs_group = factor(ttcs_group, levels = c("short","long"))),
    by = "ID"
  ) %>%
  mutate(
    stimulation = factor(stimulation, levels = stim_keep),
    marker = as.factor(marker),
    value_log = log10(pmax(value, 1e-9))  # Stabilize variance; if raw values needed, change to value_log = value
  )

## ===== 4) Calculate Hedges’ g (long − short; short as reference) + t-test + BH adjustment =====
effsize_df <- dat %>%
  group_by(stimulation, marker) %>%
  rstatix::cohens_d(value_log ~ ttcs_group,
                    var.equal = FALSE,
                    hedges.correction = TRUE,
                    ci = TRUE) %>%
  ungroup() %>%
  left_join(
    dat %>%
      group_by(stimulation, marker) %>%
      rstatix::t_test(value_log ~ ttcs_group, var.equal = FALSE) %>%
      ungroup() %>%
      mutate(p_BH = p.adjust(p, method = "BH")) %>%
      select(stimulation, marker, p_BH, statistic),
    by = c("stimulation","marker")
  ) %>%
  mutate(
    dir = case_when(statistic > 0 ~ "UP",
                    statistic < 0 ~ "DOWN",
                    TRUE ~ NA_character_),
    p_BH.signif = case_when(
      is.na(p_BH) ~ "ns",
      p_BH < 0.0001 ~ "****",
      p_BH < 0.001  ~ "***",
      p_BH < 0.01   ~ "**",
      p_BH < 0.05   ~ "*",
      TRUE ~ "ns"
    ),
    signif_dir = ifelse(p_BH.signif == "ns" | is.na(dir), "ns",
                        paste(p_BH.signif, dir))
  ) %>%
  dplyr::rename(g = effsize)  # Explicitly use dplyr::rename

## ===== 5) Sorting (by |g| global descending) & striped background rows =====
marker_order <- effsize_df %>%
  arrange(desc(abs(g))) %>%
  pull(marker) %>%
  unique() %>%
  as.character()

effsize_df <- effsize_df %>%
  mutate(marker = factor(marker, levels = marker_order))

strip_bg <- effsize_df %>%
  group_by(stimulation) %>%
  distinct(marker) %>%
  arrange(stimulation, marker) %>%
  mutate(row_id = row_number(),
         ymin   = row_id - 0.5,
         ymax   = row_id + 0.5) %>%
  filter(row_id %% 2 == 1) %>%
  ungroup()

## ===== 6) Color mapping =====
color_mapping <- c(
  "ns"       = "#888888",
  "* UP"     = "#ffd4c4",
  "** UP"    = "#FF9F94",
  "*** UP"   = "#ff5c33",
  "**** UP"  = "#cc2900",
  "* DOWN"   = "#c2d1ff",
  "** DOWN"  = "#99b3ff",
  "*** DOWN" = "#668cff",
  "**** DOWN"= "blue"
)

## ===== 7) Forest/point plot (Hedges’ g, long − short) =====
p_g <- ggplot(effsize_df, aes(x = g, y = marker)) +
  geom_rect(data = strip_bg, inherit.aes = FALSE,
            aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
            fill = "gray95", colour = NA) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.4) +
  geom_segment(aes(x = conf.low, xend = conf.high, yend = marker), size = 0.4) +
  geom_point(aes(fill = signif_dir), shape = 22, size = 3) +
  geom_text(aes(x = max(conf.high, na.rm = TRUE) + 0.1,
                label = ifelse(p_BH.signif == "ns", "", p_BH.signif)),
            hjust = 0, size = 3.5, family = "sans") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.30))) +
  scale_y_discrete(drop = TRUE, expand = expansion(add = c(0.5, 0.5))) +
  scale_fill_manual(values = color_mapping, na.value = "#888888") +
  facet_grid(. ~ stimulation, scales = "free_x") +
  labs(x = "Hedges' g (long − short), 95% CI", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(face = "bold"),
    strip.text  = element_text(face = "bold"),
    panel.grid  = element_blank(),
    axis.line.x = element_line(size = 0.25),
    legend.position = "none",
    plot.margin = margin(10, 20, 10, 10)
  )

print(p_g)
# ggsave("hedges_forest_ttcs_by_stim.svg", p_g, width = 12, height = 6, dpi = 300)


library(tidyverse)
library(rstatix)
library(ggpubr)

stim_keep <- c("KP","PP","LPS")

# Data preparation: merge groups, log10 transform
dat_log <- luminex_long %>%
  filter(stimulation %in% stim_keep) %>%
  inner_join(
    merged_data %>%
      transmute(ID = as.character(EB_id),
                ttcs_group = factor(ttcs_group, levels = c("short","long"))),
    by = "ID"
  ) %>%
  mutate(
    stimulation = factor(stimulation, levels = stim_keep),
    marker = factor(marker),
    value_log10 = log10(pmax(value, 1e-9))
  )

# Statistical tests
stat_test <- dat_log %>%
  group_by(stimulation, marker) %>%
  wilcox_test(value_log10 ~ ttcs_group) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  ungroup()

# Determine y.position for significance stars (per marker × stimulation)
stat_test <- stat_test %>%
  group_by(stimulation, marker) %>%
  mutate(y.position = max(dat_log$value_log10[dat_log$stimulation == stimulation &
                                                dat_log$marker == marker], na.rm = TRUE) + 0.3) %>%
  ungroup()

# Violin plot with significance annotation
p_violin <- ggplot(dat_log, aes(x = ttcs_group, y = value_log10, fill = ttcs_group)) +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.6) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.12, alpha = 0.3, size = 0.7) +
  facet_grid(stimulation ~ marker, scales = "free_y") +
  stat_pvalue_manual(
    stat_test,
    label = "p.adj.signif",   # Significance stars
    tip.length = 0.01,
    hide.ns = TRUE            # Do not display "ns"
  ) +
  labs(x = "TTCS group", y = "log10(Value)",
       title = "Biomarker distributions by TTCS group (log10)") +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none"
  )

print(p_violin)


#take as continous
### linear regression with "marker level after stinulation" as outcome and days to stability as predictor
#not based on logged
dat2 <- luminex_long %>%
  filter(stimulation %in% stim_keep) %>%
  inner_join(
    merged_data %>%
      transmute(ID = as.character(EB_id),
                ttcs = ttcs_halm_372_days ),
    by = "ID"
  ) %>%
  mutate(
    stimulation = factor(stimulation, levels = stim_keep),
    marker = as.factor(marker),
    value_log = log10(pmax(value, 1e-9))  # Stabilize variance; if raw values needed, change to value_log = value
  )


table(dat2$group)

dat2$stim_mark <- paste0(dat2$stimulation, "_", dat2$marker)

library(dplyr)
library(broom)


dat2 <- dat2[dat2$ttcs!=0,]

results <- dat2 %>%
  group_by(stim_mark) %>%
  do(tidy(lm(value_log ~ ttcs, data = .))) %>%
  ungroup() %>%
  filter(term == "ttcs") %>%              # keep only the slope term
  mutate(
    logFC = estimate,                     # slope (effect size)
    negLogP = -log10(p.value)             # for volcano plot
  )




#####
library(dplyr)
library(ggplot2)
library(ggrepel)

res <- dat2 %>%
  group_by(stim_mark, stimulation, marker) %>%
  do(broom::tidy(lm(value_log ~ ttcs, data = .))) %>%
  ungroup() %>%
  filter(term == "ttcs") %>%
  mutate(
    beta    = estimate,
    negLogP = -log10(p.value),
    sig     = p.value < 0.05   
  )

ggplot(res, aes(x = beta, y = negLogP, color = stimulation)) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point(alpha = 0.9) +
  geom_text_repel(
    data = res %>% filter(sig) %>% slice_max(order_by = negLogP, n = 10),
    aes(label = stim_mark),   
    size = 3, max.overlaps = 100
  ) +
  labs(x = "Change in cytokine response per day to stabilityS",
       y = "-log10(p-value)",
       title = "Volcano of TTCS associations (unadjusted p)") +
  theme_classic(base_size = 12)


##using aonova 
models <- dat2 %>%
  mutate(logttcs_fac = as.factor(logttcs)) %>%
  group_by(stim_mark) %>%
  do(model = aov(value_log ~ logttcs_fac, data = .))

# stim_mark
summary(models$model[[1]])
