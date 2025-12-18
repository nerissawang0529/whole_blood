rm(list = ls())

# ---- Packages ----
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, rio, lubridate, stringr, purrr)

# 1.---- Data ----
studydata <- import("Original_data/studydata_EB.csv") %>%
  filter(group == "CAP")


# 2.---- Long -> Wide microbiology table (tests 1–7) ----
micro_long <- studydata %>%
  select(
    EB_id, admission_dt1,
    matches("^outcome_pat_microbiology_results_1_Microbiological test [1-7]_(Sample origin|Sample date|Result|VIRUS found|BACTERIA found|FUNGI found|Specify other pathogen)$")
  ) %>%
  pivot_longer(
    cols = -c(EB_id, admission_dt1),
    names_to = c("test_no", "field"),
    names_pattern = "^outcome_pat_microbiology_results_1_Microbiological test ([1-7])_(.+)$",
    values_to = "value"
  ) %>%
  pivot_wider(names_from = field, values_from = value)


# 3.---- Keep relevant sample origins & drop obvious contamination ----
micro_filtered <- micro_long %>%
  filter(
    is.na(`Sample origin`) | str_detect(str_to_lower(`Sample origin`), "blood|sputum|bronchoalveolar"),
    is.na(Result) | Result != "Atypical result (contamination)"
  )


# 4.---- Normalize placeholders and fill from 'Specify other pathogen' when useful ----
bact_placeholders  <- c(NA, "", "Other", "No bacterial pathogen found", "Not applicable", "No growth/pathogen found")
fungi_placeholders <- c(NA, "", "No fungal pathogen found", "Not applicable")

micro_updated <- micro_filtered %>%
  mutate(
    .spec_raw = trimws(`Specify other pathogen`),
    .spec_lc  = str_to_lower(.spec_raw),
    # fungi keywords (extend as needed)
    .is_fungi = str_detect(.spec_lc, "\\baspergillus\\b"),
    # bacteria keywords (allow non-letters between genus/species)
    .is_bact  = str_detect(
      .spec_lc,
      "\\bmoraxella\\W*catarrhalis\\b|\\bserratia\\W*marcescens\\b|\\bstenotrophomonas\\W*maltophilia\\b|\\bpseudomonas\\W*fluorescens\\b"
    ),
    `FUNGI found`    = if_else(.is_fungi & (`FUNGI found` %in% fungi_placeholders), .spec_raw, `FUNGI found`),
    `BACTERIA found` = if_else(.is_bact  & (`BACTERIA found` %in% bact_placeholders),  .spec_raw, `BACTERIA found`)
  ) %>%
  select(-.spec_raw, -.spec_lc, -.is_fungi, -.is_bact)


# 5.---- Date fixes & parsing ----
micro_updated <- micro_updated %>%
  mutate(
    `Sample date` = dmy(`Sample date`),
    `Sample date` = case_when(
      EB_id == 3015 & `Sample date` == dmy("21-11-2017") ~ dmy("21-12-2017"),
      EB_id == 3031 & `Sample date` == dmy("15-07-2018") ~ dmy("15-05-2018"),
      EB_id == 1167 & `Sample date` == dmy("18-11-2021") ~ dmy("18-11-2020"),
      EB_id == 7003 & `Sample date` == dmy("26-12-2019") ~ dmy("26-11-2019"),
      TRUE ~ `Sample date`
    )
  ) %>%
  mutate(
    admission_dt1 = ymd(admission_dt1),
    # robust parsing in case formats differ
    sample_dt_try1 = suppressWarnings(ymd(`Sample date`)),
    sample_dt_try2 = suppressWarnings(dmy(`Sample date`)),
    sample_dt     = coalesce(sample_dt_try1, sample_dt_try2, `Sample date`),
    days_after_adm = as.numeric(difftime(sample_dt, admission_dt1, units = "days"))
  ) %>%
  select(-sample_dt_try1, -sample_dt_try2)


# 6.---- Clean organism fields (lowercase + trim + NA for placeholders) ----
bact_ph_lc  <- str_squish(str_to_lower(c("", "No bacterial pathogen found", "Not applicable", "Other", "No growth/pathogen found")))
fungi_ph_lc <- str_squish(str_to_lower(c("", "No fungal pathogen found", "Not applicable", "Other")))

df_clean <- micro_updated %>%
  mutate(
    bact_clean  = str_squish(str_to_lower(as.character(`BACTERIA found`))),
    fungi_clean = str_squish(str_to_lower(as.character(`FUNGI found`))),
    bact_clean  = if_else(is.na(`BACTERIA found`) | bact_clean  %in% bact_ph_lc,  NA_character_, bact_clean),
    fungi_clean = if_else(is.na(`FUNGI found`)   | fungi_clean %in% fungi_ph_lc, NA_character_, fungi_clean)
  )


# 7.---- Per-EB_id baseline vs ≥48h sets, and secondary infection flag ----
micro_updated <- df_clean %>%
  group_by(EB_id) %>%
  mutate(
    # sets at baseline (day 0) and later (≥2 days ~ ≥48h)
    bact_base   = list(unique(na.omit(bact_clean [days_after_adm == 0]))),
    fungi_base  = list(unique(na.omit(fungi_clean[days_after_adm == 0]))),
    bact_later  = list(unique(na.omit(bact_clean [days_after_adm >= 2]))),
    fungi_later = list(unique(na.omit(fungi_clean[days_after_adm >= 2]))),
    
    has_base    = lengths(bact_base)  + lengths(fungi_base)  > 0,
    has_later   = lengths(bact_later) + lengths(fungi_later) > 0,
    
    new_bact    = map2_lgl(bact_later,  bact_base,  ~ length(setdiff(.x, .y)) > 0),
    new_fungi   = map2_lgl(fungi_later, fungi_base, ~ length(setdiff(.x, .y)) > 0),
    
    # Definition: later result present AND (no baseline organism OR any novel organism appears)
    secondary_infection = has_later & (!has_base | new_bact | new_fungi)
  ) %>%
  ungroup()


# 8.---- Final compact output ----
micro_updated_2 <- micro_updated %>%
  select(
    EB_id, admission_dt1, `Sample origin`, `Sample date`,
    bact_clean, fungi_clean, bact_base, fungi_base, bact_later, fungi_later,
    secondary_infection
  )


#In total, 236 CAP patients were analyzed, of whom 1004, 1008, 1098, 1101, 3015, and 3014 had secondary infection.
