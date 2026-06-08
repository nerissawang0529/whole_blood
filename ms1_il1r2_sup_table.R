rm(list = ls())

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tableone)
  library(tibble)
  library(tools)
  library(readxl)
})

############################################################
## Global formatting options
############################################################
CONT_DIGITS <- 1
CAT_DIGITS  <- 1
P_DIGITS    <- 3

pacman::p_load(
  pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
  nlme, lcmm, tableone, lattice, reshape2, data.table, scales, Hmisc,
  kmed, survival, survminer, skimr, grDevices, ggrepel
)

############################################################
## Load data
############################################################
immunocompromised <- read_excel("Original_data/EB_corrected_immunocompromised.xlsx")

studydata <- read_csv(
  "Original_data/studydata_with_ms1_il1r2.csv",
  show_col_types = FALSE
)

############################################################
## Merge immunocompromised information into studydata
############################################################
immunocompromised_keep <- immunocompromised %>%
  mutate(
    EB_id = str_extract(as.character(ID), "^\\d+")
  ) %>%
  select(
    EB_id,
    IMM_COMPROM_JN,
    HIV_JN,
    NEUTROPENIE_JN,
    IMM_SUPPRES_JN,
    CHEMO_JN,
    TRANSPLANTATIE_JN,
    LEUKEMIE_JN
  ) %>%
  mutate(
    across(
      c(
        IMM_COMPROM_JN,
        HIV_JN,
        NEUTROPENIE_JN,
        IMM_SUPPRES_JN,
        CHEMO_JN,
        TRANSPLANTATIE_JN,
        LEUKEMIE_JN
      ),
      ~ factor(.x, levels = c("No", "Yes"))
    )
  )

studydata <- studydata %>%
  mutate(EB_id = as.character(EB_id)) %>%
  left_join(immunocompromised_keep, by = "EB_id")

############################################################
## Variables
############################################################
allvars <- c(
  "IL1R2_group", "MS1_group",
  "age_yrs", "gender", "ethnic_group#White",
  "BMI", "symptom_days",
  "CCI", "hypertension", "cpd", "COPD", "diabetes",
  "ccd", "ckd", "mneoplasm", "immune_sup", "cnd",
  "IMM_COMPROM_JN", "IMM_SUPPRES_JN", "TRANSPLANTATIE_JN",
  "LEUKEMIE_JN", "CHEMO_JN", "HIV_JN", "NEUTROPENIE_JN",
  "qSOFA_score", "MEWS_score", "CURB_score", "PSI_new",
  "WBC_2_1", "Neutrophil_unit_1", "Lymphocyte_1_1", "lab_pat_monocytes_1",
  "Creatinine_value_1", "Platelets_value_1",
  "Blood_Urea_Nitrogen_value_1",
  "length_of_stay",
  "Non_invasive_ventilation_1", "Invasive_ventilation_1",
  "lenght_of_intubation",
  "mortality_d30", "mortality_d90", "ttcs_halm_372_days"
)

catvars <- c(
  "IL1R2_group", "MS1_group", "gender",
  "ethnic_group#White", "hypertension", "cpd", "COPD", "diabetes",
  "ccd", "ckd", "mneoplasm", "immune_sup", "cnd",
  "IMM_COMPROM_JN", "IMM_SUPPRES_JN", "TRANSPLANTATIE_JN",
  "LEUKEMIE_JN", "CHEMO_JN", "HIV_JN", "NEUTROPENIE_JN",
  "Non_invasive_ventilation_1", "Invasive_ventilation_1",
  "mortality_d30", "mortality_d90"
)

nonnormal <- c(
  "age_yrs", "BMI", "symptom_days", "CCI",
  "qSOFA_score", "MEWS_score", "CURB_score", "PSI_new",
  "WBC_2_1", "Neutrophil_unit_1", "Lymphocyte_1_1", "lab_pat_monocytes_1",
  "Creatinine_value_1", "Platelets_value_1",
  "Blood_Urea_Nitrogen_value_1",
  "length_of_stay",
  "lenght_of_intubation", "ttcs_halm_372_days"
)

drop_vars <- c("Non_invasive_ventilation_1", "Invasive_ventilation_1")

allvars   <- setdiff(allvars, drop_vars)
catvars   <- setdiff(catvars, drop_vars)
nonnormal <- setdiff(nonnormal, drop_vars)

############################################################
## TableOne -> df
############################################################
tableone_to_df <- function(data, strata_var, allvars, catvars, nonnormal) {
  
  vars_use   <- setdiff(allvars, strata_var)
  factor_use <- intersect(setdiff(catvars, strata_var), vars_use)
  
  tab <- CreateTableOne(
    vars       = vars_use,
    data       = data,
    strata     = strata_var,
    factorVars = factor_use,
    test       = TRUE
  )
  
  mat <- print(
    tab,
    nonnormal   = intersect(nonnormal, vars_use),
    quote       = FALSE,
    noSpaces    = TRUE,
    smd         = TRUE,
    missing     = TRUE,
    printToggle = FALSE,
    catDigits   = CAT_DIGITS,
    contDigits  = CONT_DIGITS,
    pDigits     = P_DIGITS
  )
  
  df <- as.data.frame(mat, check.names = FALSE)
  df$Variable <- rownames(df)
  rownames(df) <- NULL
  
  df %>% relocate(Variable)
}

############################################################
## Relabel patterns
############################################################
relabel_variable <- function(x) {
  
  x0 <- x
  x  <- str_replace_all(x, "[[:space:]]*=[[:space:]]*", "=")
  xl <- tolower(x)
  
  if (x == "n") return("n")
  
  ##########################################################
  ## Top grouping rows
  ##########################################################
  if (str_detect(x, "^IL1R2_group") && str_detect(x, "IL1R2_positive")) {
    return("IL1R2 group, IL1R2 positive (%)")
  }
  
  if (str_detect(x, "^MS1_group") && str_detect(x, "MS1_high")) {
    return("MS1 group, MS1 high (%)")
  }
  
  ##########################################################
  ## Demographics
  ##########################################################
  if (str_detect(x, "^age_yrs")) return("Age, years, median [IQR]")
  
  if (str_detect(x, "^gender") && str_detect(xl, "male")) {
    return("Sex, male, n (%)")
  }
  
  if (str_detect(x, "^ethnic_group#White") && str_detect(xl, "(=1|=true|white)")) {
    return("Ethnic group, White (%)")
  }
  
  if (str_detect(x, "^BMI")) return("Body mass index, median [IQR]")
  if (str_detect(x, "^symptom_days")) return("Duration of symptoms, days (median [IQR])")
  
  ##########################################################
  ## Comorbidities
  ##########################################################
  if (str_detect(x, "^CCI")) return("Charlson comorbidity index (median [IQR])")
  if (str_detect(x, "^hypertension")) return("Hypertension (%)")
  if (str_detect(x, "^cpd")) return("Chronic pulmonary disease (%)")
  if (str_detect(x, "^COPD")) return("COPD (%)")
  if (str_detect(x, "^diabetes")) return("Diabetes (%)")
  if (str_detect(x, "^ccd")) return("Cerebrovascular disease (%)")
  if (str_detect(x, "^ckd")) return("Chronic kidney disease (%)")
  if (str_detect(x, "^mneoplasm")) return("(Prior) malignancy (%)")
  if (str_detect(x, "^immune_sup")) return("Immunosuppressed (%)")
  if (str_detect(x, "^cnd")) return("Congestive heart disease (%)")
  
  ##########################################################
  ## Immune deficiency variables
  ##########################################################
  if (str_detect(x, "^IMM_COMPROM_JN") && str_detect(x, "Yes")) {
    return("Immune deficiency")
  }
  
  if (str_detect(x, "^IMM_SUPPRES_JN") && str_detect(x, "Yes")) {
    return("Immunosuppressives")
  }
  
  if (str_detect(x, "^TRANSPLANTATIE_JN") && str_detect(x, "Yes")) {
    return("Solid organ or bone-marrow transplant")
  }
  
  if (str_detect(x, "^LEUKEMIE_JN") && str_detect(x, "Yes")) {
    return("Haematological malignancy")
  }
  
  if (str_detect(x, "^CHEMO_JN") && str_detect(x, "Yes")) {
    return("Chemotherapy <6 months")
  }
  
  if (str_detect(x, "^HIV_JN") && str_detect(x, "Yes")) {
    return("HIV infection")
  }
  
  if (str_detect(x, "^NEUTROPENIE_JN") && str_detect(x, "Yes")) {
    return("Neutropenia")
  }
  
  ##########################################################
  ## Severity
  ##########################################################
  if (str_detect(x, "^qSOFA_score")) return("qSOFA")
  if (str_detect(x, "^MEWS_score")) return("MEWS")
  if (str_detect(x, "^CURB_score")) return("CURB-65")
  if (str_detect(x, "^PSI_new")) return("Pneumonia severity index")
  
  ##########################################################
  ## Labs
  ##########################################################
  if (str_detect(x, "^WBC_2_1")) return("Leukocyte counts, ×10⁹/L")
  if (str_detect(x, "^Neutrophil_unit_1")) return("Neutrophil counts, ×10⁹/L")
  if (str_detect(x, "^Lymphocyte_1_1")) return("Lymphocyte counts, ×10⁹/L")
  if (str_detect(x, "^lab_pat_monocytes_1")) return("Monocyte counts, ×10⁹/L")
  if (str_detect(x, "^Creatinine_value_1")) return("Creatinine, μmol/L")
  if (str_detect(x, "^Platelets_value_1")) return("Platelets count, ×10⁹/L")
  if (str_detect(x, "^Blood_Urea_Nitrogen_value_1")) return("Urea, mmol/L")
  
  ##########################################################
  ## Outcome
  ##########################################################
  if (str_detect(x, "^length_of_stay")) {
    return("Length of hospital stay, day (median [IQR])")
  }
  
  if (str_detect(x, "^mortality_d30")) return("30-d mortality (%)")
  if (str_detect(x, "^mortality_d90")) return("90-d mortality (%)")
  
  if (str_detect(x, "^ttcs_halm_372_days")) {
    return("Time to clinical stability, days (median [IQR])")
  }
  
  return(x0)
}

############################################################
## Section builder
############################################################
format_with_sections <- function(df) {
  
  df <- df %>%
    mutate(Variable = vapply(Variable, relabel_variable, character(1)))
  
  top_vars <- c(
    "n",
    "IL1R2 group, IL1R2 positive (%)",
    "MS1 group, MS1 high (%)"
  )
  
  demo_vars <- c(
    "Age, years, median [IQR]",
    "Sex, male, n (%)",
    "Ethnic group, White (%)",
    "Body mass index, median [IQR]",
    "Duration of symptoms, days (median [IQR])"
  )
  
  comorb_vars <- c(
    "Charlson comorbidity index (median [IQR])",
    "Hypertension (%)",
    "Chronic pulmonary disease (%)",
    "COPD (%)",
    "Diabetes (%)",
    "Cerebrovascular disease (%)",
    "Chronic kidney disease (%)",
    "(Prior) malignancy (%)",
    "Congestive heart disease (%)",
    "Immune deficiency",
    "Immunosuppressives",
    "Solid organ or bone-marrow transplant",
    "Haematological malignancy",
    "Chemotherapy <6 months",
    "HIV infection",
    "Neutropenia"
  )
  
  severity_vars <- c(
    "qSOFA",
    "MEWS",
    "CURB-65",
    "Pneumonia severity index"
  )
  
  lab_vars <- c(
    "Leukocyte counts, ×10⁹/L",
    "Neutrophil counts, ×10⁹/L",
    "Lymphocyte counts, ×10⁹/L",
    "Monocyte counts, ×10⁹/L",
    "Creatinine, μmol/L",
    "Platelets count, ×10⁹/L",
    "Urea, mmol/L"
  )
  
  outcome_vars <- c(
    "Length of hospital stay, day (median [IQR])",
    "30-d mortality (%)",
    "90-d mortality (%)",
    "Time to clinical stability, days (median [IQR])"
  )
  
  keep_all <- c(
    top_vars,
    demo_vars,
    comorb_vars,
    severity_vars,
    lab_vars,
    outcome_vars
  )
  
  df <- df %>%
    filter(Variable %in% keep_all)
  
  ##########################################################
  ## Important fix:
  ## Some subsets have zero HIV / neutropenia cases.
  ## Then TableOne may not print the "Yes" row.
  ## This forces all expected rows to remain in every table.
  ##########################################################
  missing_vars <- setdiff(keep_all, df$Variable)
  
  if (length(missing_vars) > 0) {
    add_df <- as.data.frame(
      matrix("", nrow = length(missing_vars), ncol = ncol(df)),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    names(add_df) <- names(df)
    add_df$Variable <- missing_vars
    
    df <- bind_rows(df, add_df)
  }
  
  header_row <- function(title, coln) {
    x <- as.list(rep("", length(coln)))
    names(x) <- coln
    x$Variable <- title
    as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
  }
  
  coln <- names(df)
  
  out <- bind_rows(
    df %>%
      filter(Variable %in% top_vars) %>%
      mutate(Variable = factor(Variable, levels = top_vars)) %>%
      arrange(Variable) %>%
      mutate(Variable = as.character(Variable)),
    
    header_row("Demographics", coln),
    
    df %>%
      filter(Variable %in% demo_vars) %>%
      mutate(Variable = factor(Variable, levels = demo_vars)) %>%
      arrange(Variable) %>%
      mutate(Variable = as.character(Variable)),
    
    header_row("Comorbidities", coln),
    
    df %>%
      filter(Variable %in% comorb_vars) %>%
      mutate(Variable = factor(Variable, levels = comorb_vars)) %>%
      arrange(Variable) %>%
      mutate(Variable = as.character(Variable)),
    
    header_row("Severity on admission, median [IQR]", coln),
    
    df %>%
      filter(Variable %in% severity_vars) %>%
      mutate(Variable = factor(Variable, levels = severity_vars)) %>%
      arrange(Variable) %>%
      mutate(Variable = as.character(Variable)),
    
    header_row("Routine laboratory markers, median [IQR]", coln),
    
    df %>%
      filter(Variable %in% lab_vars) %>%
      mutate(Variable = factor(Variable, levels = lab_vars)) %>%
      arrange(Variable) %>%
      mutate(Variable = as.character(Variable)),
    
    header_row("Outcome", coln),
    
    df %>%
      filter(Variable %in% outcome_vars) %>%
      mutate(Variable = factor(Variable, levels = outcome_vars)) %>%
      arrange(Variable) %>%
      mutate(Variable = as.character(Variable))
  )
  
  out
}

############################################################
## Top rows creator
############################################################
make_top_rows <- function(data, final_cols) {
  
  fmt_np <- function(n, denom) {
    if (is.na(n) || is.na(denom) || denom == 0) return("")
    sprintf("%d (%.1f)", as.integer(n), 100 * as.numeric(n) / as.numeric(denom))
  }
  
  N_total <- nrow(data)
  
  n_ms1_low  <- sum(data$MS1_group == "MS1_low", na.rm = TRUE)
  n_ms1_high <- sum(data$MS1_group == "MS1_high", na.rm = TRUE)
  n_il1_zero <- sum(data$IL1R2_group == "IL1R2_zero", na.rm = TRUE)
  n_il1_pos  <- sum(data$IL1R2_group == "IL1R2_positive", na.rm = TRUE)
  
  n_ms1_low_il1_zero  <- sum(data$MS1_group == "MS1_low"  & data$IL1R2_group == "IL1R2_zero", na.rm = TRUE)
  n_ms1_low_il1_pos   <- sum(data$MS1_group == "MS1_low"  & data$IL1R2_group == "IL1R2_positive", na.rm = TRUE)
  n_ms1_high_il1_zero <- sum(data$MS1_group == "MS1_high" & data$IL1R2_group == "IL1R2_zero", na.rm = TRUE)
  n_ms1_high_il1_pos  <- sum(data$MS1_group == "MS1_high" & data$IL1R2_group == "IL1R2_positive", na.rm = TRUE)
  
  tab <- table(data$MS1_group, data$IL1R2_group)
  
  pval <- NA
  
  suppressWarnings({
    chi <- try(chisq.test(tab, correct = FALSE), silent = TRUE)
    
    if (!inherits(chi, "try-error") && all(chi$expected >= 5)) {
      pval <- chi$p.value
    } else {
      pval <- fisher.test(tab)$p.value
    }
  })
  
  ptxt <- format.pval(pval, digits = P_DIGITS, eps = 0.001)
  
  rn <- setNames(as.list(rep("", length(final_cols))), final_cols)
  rn$Variable <- "n"
  rn$`MS1 low group`        <- as.character(n_ms1_low)
  rn$`MS1 high group`       <- as.character(n_ms1_high)
  rn$`IL1R2 zero group`     <- as.character(n_il1_zero)
  rn$`IL1R2 positive group` <- as.character(n_il1_pos)
  
  r_il1 <- setNames(as.list(rep("", length(final_cols))), final_cols)
  r_il1$Variable <- "IL1R2 group, IL1R2 positive (%)"
  r_il1$`IL1R2 zero group`     <- fmt_np(0, N_total)
  r_il1$`IL1R2 positive group` <- fmt_np(n_il1_pos, N_total)
  r_il1$`MS1 low group`        <- fmt_np(n_ms1_low_il1_pos, n_ms1_low)
  r_il1$`MS1 high group`       <- fmt_np(n_ms1_high_il1_pos, n_ms1_high)
  r_il1$`MS1 p-value`          <- ptxt
  r_il1$`IL1R2 p-value`        <- ptxt
  
  r_ms1 <- setNames(as.list(rep("", length(final_cols))), final_cols)
  r_ms1$Variable <- "MS1 group, MS1 high (%)"
  r_ms1$`MS1 low group`        <- fmt_np(0, N_total)
  r_ms1$`MS1 high group`       <- fmt_np(n_ms1_high, N_total)
  r_ms1$`IL1R2 zero group`     <- fmt_np(n_ms1_high_il1_zero, n_il1_zero)
  r_ms1$`IL1R2 positive group` <- fmt_np(n_ms1_high_il1_pos, n_il1_pos)
  r_ms1$`MS1 p-value`          <- ptxt
  r_ms1$`IL1R2 p-value`        <- ptxt
  
  bind_rows(
    as.data.frame(rn, stringsAsFactors = FALSE, check.names = FALSE),
    as.data.frame(r_il1, stringsAsFactors = FALSE, check.names = FALSE),
    as.data.frame(r_ms1, stringsAsFactors = FALSE, check.names = FALSE)
  )
}

############################################################
## Main function
############################################################
make_dual_table1 <- function(data, prefix) {
  
  ms1_df <- tableone_to_df(data, "MS1_group", allvars, catvars, nonnormal) %>%
    rename_with(~ paste0("MS1_", .x), -Variable)
  
  il1_df <- tableone_to_df(data, "IL1R2_group", allvars, catvars, nonnormal) %>%
    rename_with(~ paste0("IL1R2_", .x), -Variable)
  
  merged <- bind_cols(ms1_df, il1_df %>% select(-Variable)) %>%
    select(
      Variable,
      matches("^MS1_") & !matches("test$|SMD$|Missing$"),
      matches("^IL1R2_") & !matches("test$|SMD$|Missing$")
    )
  
  ms1_low_col <- names(merged)[
    str_detect(names(merged), "^MS1_") &
      str_detect(names(merged), "MS1_low$")
  ]
  
  ms1_high_col <- names(merged)[
    str_detect(names(merged), "^MS1_") &
      str_detect(names(merged), "MS1_high$")
  ]
  
  il1_zero_col <- names(merged)[
    str_detect(names(merged), "^IL1R2_") &
      str_detect(names(merged), "IL1R2_zero$")
  ]
  
  il1_pos_col <- names(merged)[
    str_detect(names(merged), "^IL1R2_") &
      str_detect(names(merged), "IL1R2_positive$")
  ]
  
  final_df <- merged %>%
    rename(
      `MS1 low group`        = all_of(ms1_low_col),
      `MS1 high group`       = all_of(ms1_high_col),
      `MS1 p-value`          = MS1_p,
      `IL1R2 zero group`     = all_of(il1_zero_col),
      `IL1R2 positive group` = all_of(il1_pos_col),
      `IL1R2 p-value`        = IL1R2_p
    ) %>%
    select(
      Variable,
      `MS1 low group`, `MS1 high group`, `MS1 p-value`,
      `IL1R2 zero group`, `IL1R2 positive group`, `IL1R2 p-value`
    ) %>%
    mutate(across(everything(), ~ ifelse(is.na(.x), "", .x)))
  
  final_df <- final_df %>%
    filter(Variable != "n")
  
  top_rows <- make_top_rows(data, names(final_df))
  
  out <- bind_rows(top_rows, final_df) %>%
    mutate(across(everything(), ~ ifelse(is.na(.x), "", .x)))
  
  out <- format_with_sections(out) %>%
    mutate(across(everything(), ~ ifelse(is.na(.x), "", .x)))
  
  out_file <- paste0("Original_data/", prefix, "_table1_ms1_il1r2.csv")
  
  write.csv(
    out,
    out_file,
    row.names = FALSE,
    na = "",
    quote = TRUE
  )
  
  message("✅ Exported: ", out_file)
  
  invisible(out)
}

############################################################
## Create Table 1 for full cohort
############################################################
make_dual_table1(studydata, prefix = "main_full_cohort")

############################################################
## Create Table 1 for whole-blood stimuli subset
############################################################
wb_long <- read_csv(
  "Original_data/whole_blood_stimuli_long.csv",
  show_col_types = FALSE
)

if ("HiIL6" %in% names(wb_long)) {
  wb_long$HiIL6 <- NULL
}

stimuli_ids <- wb_long %>%
  distinct(ID) %>%
  pull(ID) %>%
  as.character()

studydata_stimuli <- studydata %>%
  mutate(EB_id = as.character(EB_id)) %>%
  filter(EB_id %in% stimuli_ids)

message("Full cohort n = ", nrow(studydata))
message("WB stimuli subset n = ", nrow(studydata_stimuli))

make_dual_table1(studydata_stimuli, prefix = "sub_wb_stimuli")

############################################################
## Create Table 1 for baseline marker subset
############################################################
bio_raw <- read.table(
  "Original_Data/EB_optimact_final_biomarkers_luminex_cba.csv",
  header = TRUE,
  sep = ";",
  dec = ",",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

bio_raw <- bio_raw %>%
  dplyr::filter(str_starts(id, "EB_")) %>%
  mutate(ID = str_remove(id, "^EB_"))

baseline_ids <- bio_raw %>%
  distinct(ID) %>%
  pull(ID) %>%
  as.character()

studydata_baseline <- studydata %>%
  mutate(EB_id = as.character(EB_id)) %>%
  filter(EB_id %in% baseline_ids)

message("Baseline marker subset n = ", nrow(studydata_baseline))

make_dual_table1(studydata_baseline, prefix = "sub_baseline_markers")

############################################################
## Create Table 1 for MS1data_complete subset
############################################################
ms1data_complete <- read_csv(
  "Original_data/MS1data_complete.csv",
  show_col_types = FALSE
)

if (!"file" %in% names(ms1data_complete)) {
  stop("The file MS1data_complete.csv must contain a column named 'file'.")
}

ms1_complete_ids <- ms1data_complete %>%
  mutate(ID = str_extract(file, "(?<=Norm_)\\d+")) %>%
  filter(!is.na(ID)) %>%
  distinct(ID) %>%
  pull(ID) %>%
  as.character()

studydata_ms1complete <- studydata %>%
  mutate(EB_id = as.character(EB_id)) %>%
  filter(EB_id %in% ms1_complete_ids)

message("MS1data_complete unique IDs found = ", length(ms1_complete_ids))
message("MS1data_complete subset n = ", nrow(studydata_ms1complete))

make_dual_table1(studydata_ms1complete, prefix = "sub_MS1data_complete")