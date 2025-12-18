rm(list = ls())

## Healthy volunteers Nerissa U219 / HTA data 

# -----------------------------------------------------------------------------#
# ------------ Compare healthies for coconut ----------------------------------#
# -----------------------------------------------------------------------------#

# Load packages 
library(dplyr)
library(tableone)
library(lubridate)
library(DESeq2)
library(rio)
library(GEOquery)
library(ggplot2)

# U219 ------------------------------------------------------------------------
load("Original_data/MARS_U219_correct_age.eset")

# inspection of healthy volunteers=> 42 available 
U219_hv <- pData(U219_correct_ages) %>%
  dplyr::filter(health == "healthy") %>%
  select(MARSID, Patient_gender, age)

# filter out healthies from entire eset 
healthy_eset_U219 <- U219_correct_ages[, U219_correct_ages$health == "healthy"]
dim(healthy_eset_U219) # correct, 42 samples 



















































# HTA ------------------------------------------------------------------------

# load series and platform data from GEO (primary set Brendon; takes a while to load)
gset <- GEOquery::getGEO("GSE134364", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]

# Inspect HV 
HTA_hv <- pData(gset) %>%
  dplyr::filter(characteristics_ch1.2 == "disease state: healthy") %>%
  dplyr::select(title, `gender:ch1`, `age:ch1`) %>%
  dplyr::rename(id = title, sex = `gender:ch1`, age = `age:ch1`) %>%
  dplyr::mutate(
    platform = "HTA",
    age = suppressWarnings(as.numeric(age)),
    sex = tolower(trimws(sex)),
    sex = dplyr::case_when(
      sex %in% c("male", "m") ~ "male",
      sex %in% c("female", "f") ~ "female",
      TRUE ~ NA_character_
    )
  )

# filter out healthies from entire eset 
healthy_eset_HTA <- gset[, gset$characteristics_ch1.2 == "disease state: healthy"]
dim(healthy_eset_HTA) # correct, 83 samples 

