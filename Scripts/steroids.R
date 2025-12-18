#the definition:
#1.new pathogen at 48 hours, no pathogen at the beginning
#2.new pathogen at 48 hours, different pathogen at the beginning
#not in the urine



rm(list = ls ())

# Make table one
## install.packages("pacman")
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
               tableone, kmed, survival, survminer, skimr, lattice, reshape2, grDevices, ggrepel)
studydata <- rio::import("Original_data/studydata.csv")
dexa_df <- studydata[, c(
  "EB_id",
  "admission_dt1",
  "base_sample_pat_blood",
  "dexa",
  "dexa_dose",
  "dexa_freq",
  "dexa_route",
  "dexa_start",
  "dexa_stop",
  "dexa_moment",
  "dexa_time",
  "group"
)]

library(dplyr)

library(dplyr)

summary_table <- dexa_df %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(
    n_total = dplyr::n(),                                 # 总样本数
    n_dexa_yes = sum(dexa == "Yes", na.rm = TRUE)         # dexa == Yes 的数量
  ) %>%
  dplyr::mutate(
    pct_dexa_yes = n_dexa_yes / n_total * 100             # 百分比
  )

print(summary_table)
