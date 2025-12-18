rm(list = ls ())

# Make table one
## install.packages("pacman")
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
               tableone, kmed, survival, survminer, skimr, lattice, reshape2, grDevices, ggrepel)
studydata <- import("Original_data/ELDER-BIOME_excel_export_20230418012733.xlsx")
## p_load(tableone)

##
list_of_patients <- readxl::read_xlsx("Original_data/2310_Present_WB_stimulation_list_EB_V3_1101_from_Stijn.xlsx")
table(list_of_patients$Stimulations)


## Stijn found new samples 2065, 3146, 1076, the experiment made a mistake 2098, "1090" "3004" "3137" are not in the luminex
studydata <- studydata[
  (studydata$`Participant Id` %in% list_of_patients$Patient_study_ID | 
     studydata$`Participant Id` %in% c(2065, 3146, 1076)) &
    !(studydata$`Participant Id` %in% c(2098, 1090, 3004, 3137)), 
]
studydata <- rbind(studydata, c(2125, rep(NA, ncol(studydata) - 1)))#2125 is not in the studydata, so rbind it in studydata

#remove the ones not retesed
ids_to_remove <- c(3024, 1194, 1195, 1210, 1227, 1228, 3095, 3098, 3099, 3104, 
                   3115, 3128, 3165, 3166, 3167, 3181, 3186, 3201, 3207, 3208, 
                   3218, 3219, 3227, 3245, 3286, 1149, 1160, 1186, 3174, 3291, 3248)

studydata <- studydata %>%
  filter(!`Participant Id` %in% ids_to_remove)


## check <- list_of_patients[!list_of_patients$`Patient study ID` %in% data$`Participant Id`, ]
lunimex_not_in_studydata <- setdiff(luminex_updated_3$ID, studydata$`Participant Id`)
studydata_not_in_lunimex <- setdiff(studydata$`Participant Id`, luminex_updated_3$ID)
#luminex_updated_3 is from code luminex_cleaned

colnames(studydata)[colnames(studydata) == "Participant Id"] <- "Record Id"


## Fix missing
studydata[studydata == "##USER_MISSING_99##"] <- NA
studydata[studydata == "##USER_MISSING_97##"] <- NA
studydata[studydata == "##USER_MISSING_96##"] <- NA
studydata[studydata == "##USER_MISSING_95##"] <- NA
studydata[studydata == "##USER_MISSING_98##"] <- NA
studydata[studydata == ""] <- NA
studydata[studydata == "Unknown"] <- NA
studydata[studydata == -99] <- NA
studydata[studydata == -97] <- NA
studydata[studydata == -96] <- NA
studydata$outcome_pat_icu[studydata$outcome_pat_icu == ""] <- NA
studydata$outcome_pat_invas_vent[studydata$outcome_pat_invas_vent == ""] <- NA
studydata$outcome_pat_invas_vent[studydata$outcome_pat_invas_vent == "98"] <- NA
studydata$outcome_pat_icu[studydata$outcome_pat_icu == "98"] <- NA


## admission date as date
studydata$base_date_admission <- as.Date(studydata$base_date_admission, format = "%d-%m-%Y")
## studydata$`Record Id` <- as.numeric(studydata$`Record Id`)

## studydata
##
## find people who are samples on the ICU ## 32 patients
studydata$outcome_pat_icu_admission <- as.Date(studydata$outcome_pat_icu_admission, "%d-%m-%Y")
class(studydata$outcome_pat_icu_admission)

studydata$base_date_assessment <- as.Date(studydata$base_date_assessment, "%d-%m-%Y")
class(studydata$base_date_assessment)

studydata$time_sampling_ICU_admission <- difftime(studydata$outcome_pat_icu_admission, studydata$base_date_assessment, units = "days")
class(studydata$time_sampling_ICU_admission)

studydata$time_ICU_admission <- difftime(studydata$outcome_pat_icu_admission, studydata$base_date_admission, units = "days")
studydata$time_ICU_admission <- as.numeric(studydata$time_ICU_admission)
studydata$sampling_on_ICU_day <- ifelse(studydata$base_date_assessment > studydata$outcome_pat_icu_admission, "Yes",
                                        ifelse(studydata$base_date_assessment == studydata$outcome_pat_icu_admission, "maybe", "No"))
table(studydata$sampling_on_ICU_day)

## manually inspect patients with ICU admission and sampling on same day
ICU_and_sampling_same_day <- filter(studydata, studydata$sampling_on_ICU_day == "maybe")
ICU_and_sampling_same_day <- select(ICU_and_sampling_same_day, "Record Id", "base_date_assessment", "outcome_pat_icu_admission", "time_sampling_ICU_admission")


studydata$sampling_on_ICU_day[studydata$`Record Id` == 1127] <- "No"
studydata$sampling_on_ICU_day[studydata$`Record Id` == 	1150] <- "No"
studydata$sampling_on_ICU_day[studydata$`Record Id` == 	1154] <- "Yes"
studydata$sampling_on_ICU_day[studydata$`Record Id` == 	1157] <- "Yes"
studydata$sampling_on_ICU_day[studydata$`Record Id` == 	1169] <- "Yes"
studydata$sampling_on_ICU_day[studydata$`Record Id` == 	1183] <- "Yes" 
studydata$sampling_on_ICU_day[studydata$`Record Id` == 	1197] <- "Yes"
studydata$sampling_on_ICU_day[studydata$`Record Id` == 	3081] <- "No"
studydata$sampling_on_ICU_day[studydata$`Record Id` == 	3147] <- "No"
studydata$sampling_on_ICU_day[studydata$`Record Id` == 	3183] <- "No"
studydata$sampling_on_ICU_day[studydata$`Record Id` == 	3189] <- "Yes"
studydata$sampling_on_ICU_day[studydata$`Record Id` == 	1209] <- "No"
studydata$sampling_on_ICU_day[studydata$`Record Id` == 	3197] <- "No"
studydata$sampling_on_ICU_day[studydata$`Record Id` == 	1227] <- "Yes"
studydata$sampling_on_ICU_day[studydata$`Record Id` == 	1091] <- "No"

## check => there should be 0 "maybe" Now
ICU_and_sampling_same_day <- filter(studydata, studydata$sampling_on_ICU_day == "maybe")
remove(ICU_and_sampling_same_day)

## make na No
studydata$ICU_patient <- ifelse(is.na(studydata$sampling_on_ICU_day), "No", studydata$sampling_on_ICU_day)
table(studydata$ICU_patient)

## make MEWS variable
studydata$base_rr <- as.numeric(studydata$base_rr)
studydata$base_oxygen <- as.numeric(studydata$base_oxygen)
studydata$base_spo2 <- as.numeric(studydata$base_spo2)
studydata$base_hr <- as.numeric(studydata$base_hr)
studydata$base_syst_bp <- as.numeric(studydata$base_syst_bp)
studydata$base_temp <- as.numeric(studydata$base_temp)
studydata$MEWS_resp_rate <- with(studydata, ifelse(is.na(base_rr), NA,
                                                   ifelse(base_rr < 9, 2,
                                                          ifelse(base_rr %in% 9:14, 0,
                                                                 ifelse(base_rr %in% 15:20, 1,
                                                                        ifelse(base_rr %in% 21:29, 2,
                                                                               ifelse(base_rr >= 30, 3, NA)))))))
###heart rate for MEWS score
studydata$MEWS_hrt_rate <- with(studydata, ifelse(is.na(base_hr), NA,
                                                  ifelse(base_hr < 40, 2,
                                                         ifelse(base_hr %in% 40:50, 1, 
                                                                ifelse(base_hr %in% 51:100, 0,
                                                                       ifelse(base_hr %in% 101:110, 1,
                                                                              ifelse(base_hr %in% 111:129, 2,
                                                                                     ifelse(base_hr >= 130, 3, NA))))))))

###systolic blood pressure for MEWS score
studydata$MEWS_sys_bp <- with(studydata, ifelse(is.na(base_syst_bp), NA,
                                                ifelse(base_syst_bp <= 70, 3,
                                                       ifelse(base_syst_bp %in% 71:80, 2, 
                                                              ifelse(base_syst_bp %in% 81:100, 1,
                                                                     ifelse(base_syst_bp %in% 101:199, 0,
                                                                            ifelse(base_syst_bp >= 200, 2, NA)))))))

###temperature for MEWS score
studydata$MEWS_temp <- with(studydata, ifelse(is.na(base_temp), NA,
                                              ifelse(base_temp < 35, 2, 
                                                     ifelse(base_temp >= 35 & base_temp < 38.5, 0,
                                                            ifelse(base_temp >= 38.5, 2, NA)))))

###level of consciousness for MEWS score
studydata$base_mental_status 
studydata$MEWS_conciousness <- ifelse(!is.na(studydata$base_mental_status) & studydata$base_mental_status == "Yes", 1, 0)


studydata$MEWS_score <- with(studydata, MEWS_resp_rate + MEWS_hrt_rate + MEWS_sys_bp +
                               MEWS_conciousness + MEWS_temp)


## CURB with missing of one variable is missing
studydata$CURB_rr  <- with(studydata, ifelse(is.na(base_rr), NA,
                                             ifelse(base_rr >= 30, 1, 0)))
studydata$CURB_confusion  <-  with(studydata, ifelse(!is.na(base_mental_status) & base_mental_status == "Yes", 1, 0))
studydata$CURB_bp  <- with(studydata, ifelse(is.na(base_syst_bp), NA,
                                             ifelse(base_syst_bp <= 90, 1, 
                                                    ifelse(base_diast <= 60, 1, 0))))
studydata$CURB_age <- with(studydata, ifelse(dem_pat_age >= 65, 1, 0))
studydata$CURB_ureum <- with(studydata, ifelse(lab_pat_bun == "No", 0,
                                               ifelse(is.na(lab_pat_bun_1), 0, 
                                                      ifelse(lab_pat_bun_1 > 7, 1, 0))))
studydata$CURB_score<- with(studydata, CURB_rr  + CURB_confusion  + CURB_bp  + CURB_age  + CURB_ureum )
studydata$CURB_score_no_age <-with(studydata, CURB_rr  + CURB_confusion  + CURB_bp  + CURB_ureum )

studydata$base_gcs <- ifelse(studydata$base_mental_status == "Yes", studydata$base_gcs, 15)

## calculate qSOFA
studydata$qSOFA_bp <- with(studydata, ifelse(is.na(base_syst_bp), NA,
                                             ifelse(base_syst_bp <= 100, 1, 0)))
studydata$qSOFA_altered_mentation <- with(studydata, ifelse(studydata$base_mental_status == "Yes" | studydata$base_gcs <15, 1, 0))
studydata$qSOFA_rr <- with(studydata, ifelse(is.na(base_rr), NA,
                                             ifelse(base_rr >= 22, 1, 0)))
studydata$qSOFA_score <- with(studydata,qSOFA_bp + qSOFA_altered_mentation + qSOFA_rr)

## make 'white race' variable
studydata$white_race <- ifelse(studydata$dem_pat_race == "Caucasian", 1, 0)

## calculate time between admission and symptoms
studydata$base_date_onset <- as.Date(studydata$base_date_onset, format = "%d-%m-%Y")
studydata$days_symptoms <- difftime(studydata$base_date_admission,studydata$base_date_onset, units = "days")
studydata$days_symptoms <- ifelse(studydata$days_symptoms <0 | studydata$days_symptoms >30, NA, studydata$days_symptoms)

## calculate time between admission and sampling
studydata$base_date_assessment
studydata$base_date_assessment <- as.Date(studydata$base_date_assessment, "%d-%m-%Y")
studydata$days_to_sampling <- difftime(studydata$base_date_assessment, studydata$base_date_admission, units = "days")

studydata$statin <- studydata$`chrmed_pat_other.Statins_simvastatin_atorvastatin`
studydata$PPI <- studydata$`chrmed_pat_other.Proton_pump_inhibitors_omeprazol_pantoprazol`
studydata$base_oxygen <- ifelse(!is.na(studydata$base_oxygen) & studydata$base_oxygen >= 1, "Yes", "No")

## create mortality like predict
table(studydata$day28_pat_status)
studydata$day28_mortality <- ifelse(studydata$day28_pat_status == "Alive, discharged to home (or previous residency or rehabilitation clinic)" 
                                    | studydata$day28_pat_status == "Alive, transferred to other hospital" | studydata$day28_pat_status == "Alive, readmitted to the hospital after earlier discharge", 0,
                                    ifelse(studydata$day28_pat_status == "Alive, still hospitalized on day 28 after admission", 0, 
                                           ifelse(studydata$day28_pat_status == "Deceased", 1, NA)))

studydata$wk12_mortality <- ifelse(studydata$day90_pat_status == "Alive, discharged to home (or previous residency or rehabilitation clinic)" 
                                   | studydata$day90_pat_status == "Alive, transferred to other hospital" | studydata$day90_pat_status == "Alive, readmitted to the hospital after earlier discharge", 0,
                                   ifelse(studydata$day90_pat_status == "Alive, still hospitalized on day 90 after admission", 0, 
                                          ifelse(studydata$day28_mortality == 1, 1,
                                                 ifelse(studydata$day90_pat_status == "Deceased", 1, NA))))

## remove strange dates
studydata$day28_pat_death[studydata$day28_pat_death == "01-01-2999"] <- NA
studydata$day28_pat_death[studydata$day28_pat_death == "01-01-2996"] <- NA
studydata$day90_pat_death[studydata$day90_pat_death == "01-01-2999"] <- NA
studydata$day90_pat_death[studydata$day90_pat_death == "01-01-2996"] <- NA

## in hospital mortality
studydata$date_of_death <- as.Date(studydata$day28_pat_death, format = "%d-%m-%Y")
studydata$day28_pat_discharge <- as.Date(studydata$day28_pat_discharge, format = "%d-%m-%Y")
studydata$hospdeath <- ifelse(studydata$day28_pat_discharge <= studydata$date_of_death, 1, 0)
studydata$hospdeath <- ifelse(is.na(studydata$hospdeath), 0, studydata$hospdeath)
table(studydata$hospdeath)

##
studydata$diabetes <- ifelse(studydata$`comorb_pat_def_metabolic#Diabetes type 1` == 1, 1, 
                             ifelse(studydata$`comorb_pat_def_metabolic#Diabetes type 2 (insulin dependent)` == 1, 1, 
                                    ifelse(studydata$`comorb_pat_def_metabolic#Diabetes type 2 (only oral medication)` == 1, 1, 0)))

## rename to biobank variables
## demographics, waves is calculated
names(studydata)[names(studydata) == "Record Id"] <- "EB_id"
names(studydata)[names(studydata) == "base_date_admission"] <- "admission_dt1"
names(studydata)[names(studydata) == "dem_pat_gender"] <- "gender"
names(studydata)[names(studydata) == "white_race"] <- "ethnic_group#White"
names(studydata)[names(studydata) == "dem_pat_age"] <- "age_yrs"
names(studydata)[names(studydata) == "dem_pat_bmi"] <- "BMI"
names(studydata)[names(studydata) == "base_date_onset"] <- "onset_dt"
names(studydata)[names(studydata) == "days_symptoms"] <- "symptoms_days"

## comorbidities => missing petpic, reumathologic = connective, hematologic cancer not split out
sum(studydata$`comorb_pat_def_liver.1` == T) ## cirrhosis
sum(studydata$`comorb_pat_def_liver.2` == T) ## hep B
sum(studydata$`comorb_pat_def_liver.3` == T) ## hep c
sum(studydata$`comorb_pat_def_liver.4` == T) ## Primary sclerosing cholangitis 
studydata$live_disease <- "No"

names(studydata)[names(studydata) == "comorb_pat_def_cvd.Past_myocardial_infarction"] <- "vasc_cvd.Myocardial_infarction_or_acute_coronary_syndrome"
names(studydata)[names(studydata) == "comorb_pat_def_cvd.Congestive_heart_failure"] <- "vasc_cvd.Chronic_heart_failure"
names(studydata)[names(studydata) == "comorb_pat_def_cvd.Peripheral_ischaemicatherosclerotic_disease"] <- "perif_vasc"
names(studydata)[names(studydata) == "comorb_pat_def_cvd.Hypertension_medicated"] <- "hypertension"
names(studydata)[names(studydata) == "comorb_pat_def_cvd.Past_cerebrovascular_disease_thrombosis_embolism"] <- "vasc_cvd.Ischemic_stroke_or_TIA"
names(studydata)[names(studydata) == "comorb_pat_def_cpd.Chronic_Obstructive_Pulmonary_Disease_COPD"] <- "COPD"
names(studydata)[names(studydata) == "comorb_pat_def_other.Rheumatologicautoimmune_disease_eg_RA_SLE_PMR"] <- "connective"
studydata$peptic <- "No"
studydata$liver_charl <- "none"
studydata$hemiplegia <- "No"
names(studydata)[names(studydata) == "comorb_pat_cpd"] <- "cpd"
names(studydata)[names(studydata) == "comorb_pat_cvd"] <- "ccd"
names(studydata)[names(studydata) == "comorb_pat_def_other.Chronic_renal_disease"] <- "ckd"
names(studydata)[names(studydata) == "comorb_pat_malignancy"] <- "mneoplasm" 
names(studydata)[names(studydata) == "comorb_pat_tumor_metastatic"] <- "metastase"
names(studydata)[names(studydata) == "comorb_pat_neuro"] <- "cnd"
names(studydata)[names(studydata) == "comorb_pat_immunosup"] <- "immune_sup"
studydata$COPD <- ifelse(studydata$`comorb_pat_def_cpd#Chronic Obstructive Pulmonary Disease (COPD)` == 1, "Yes", "No")
studydata$connective <- ifelse(studydata$`comorb_pat_def_other#Rheumatologic/autoimmune disease (e.g. RA, SLE, PMR)` == 1, "Yes", "No")
studydata$perif_vasc <- ifelse(studydata$`comorb_pat_def_cvd#Peripheral ischaemic/atherosclerotic disease` == 1, "Yes", "No")

## vital signs
names(studydata)[names(studydata) == "base_temp"] <- "Temperature"
names(studydata)[names(studydata) == "base_hr"] <- "HtR"
names(studydata)[names(studydata) == "base_syst_bp"] <- "sys_bp"
names(studydata)[names(studydata) == "base_diast"] <- "dias_bp"
names(studydata)[names(studydata) == "base_rr"] <- "rtr"
names(studydata)[names(studydata) == "base_spo2"] <- "oxygen_saturation"
names(studydata)[names(studydata) == "base_gcs"] <- "EMV"

## severity scores, 4C, qsofa, mews, CURB is calculated
names(studydata)[names(studydata) == "base_pat_PSI"] <- "PSI_score"

## lab
names(studydata)[names(studydata) == "lab_pat_def_hemoglobin_1"] <- "Haemoglobin_value_1"

names(studydata)[names(studydata) == "lab_pat_def_crp"] <- "crp_1_1"
names(studydata)[names(studydata) == "lab_pat_procalcitonin_1"] <- "procalcitonine"
names(studydata)[names(studydata) == "lab_pat_lactate_1"] <- "Lactate_2_1"
names(studydata)[names(studydata) == "lab_pat_ldh_1"] <- "LDH"

names(studydata)[names(studydata) == "lab_pat_wbc_1"] <- "WBC_2_1"
names(studydata)[names(studydata) == "lab_pat_neutrophils_1"] <- "Neutrophil_unit_1"
names(studydata)[names(studydata) == "lab_pat_lymphocytes_1"] <- "Lymphocyte_1_1"

names(studydata)[names(studydata) == "lab_pat_creatinine_1"] <- "Creatinine_value_1"
names(studydata)[names(studydata) == "lab_pat_bun_1"] <- "Blood_Urea_Nitrogen_value_1"
names(studydata)[names(studydata) == "lab_pat_sodium_1"] <- "Sodium_1_1"
names(studydata)[names(studydata) == "lab_pat_potassium_1"] <- "Potassium_1_1"
names(studydata)[names(studydata) == "Glucose_admission_value"] <- "Glucose_unit_1_1"

names(studydata)[names(studydata) == "lab_pat_asat_1"] <- "AST_SGOT_1_1"
names(studydata)[names(studydata) == "lab_pat_alat_1"] <- "ALT_SGPT_1_1"
names(studydata)[names(studydata) == "lab_pat_bilirubin_1"] <- "Total_Bilirubin_2_1"

names(studydata)[names(studydata) == "lab_pat_def_trombos"] <- "Platelets_value_1"
names(studydata)[names(studydata) == "lab_pat_ddimer_1"] <- "d_dimer"
names(studydata)[names(studydata) == "lab_pat_pt_1"] <- "pt_spec"
names(studydata)[names(studydata) == "lab_pat_aptt_1"] <- "APT_APTR_1_1"

names(studydata)[names(studydata) == "lab_pat_bloodgas_pH"] <- "PH_value_1"
names(studydata)[names(studydata) == "lab_pat_bloodgas_paO2"] <- "PaO2_1"
names(studydata)[names(studydata) == "lab_pat_bloodgas_pCO2"] <- "PCO2_1"

names(studydata)[names(studydata) == "outcome_pat_supo2"] <- "Oxygen_therapy_1"
names(studydata)[names(studydata) == "outcome_pat_supo2_start"] <- "oxygen_start"
names(studydata)[names(studydata) == "outcome_pat_supo2_stop"] <- "oxygen_stop"

names(studydata)[names(studydata) == "day28_pat_death"] <- "death_date"
names(studydata)[names(studydata) == "day28_pat_discharge"] <- "discharge_date_local"

names(studydata)[names(studydata) == "outcome_pat_def_complication.Pulmonary_embolism"] <- "LE"
names(studydata)[names(studydata) == "outcome_pat_def_complication.Peripheral_venous_thrombosis"] <- "DVT"
names(studydata)[names(studydata) == "outcome_pat_def_complication.Bacteremia"] <- "bacteremia"
names(studydata)[names(studydata) == "outcome_pat_def_complication.Acute_lung_injury__ARDS"] <- "ARDS"



names(studydata)[names(studydata) == "outcome_pat_icu"] <- "ICU_Medium_Care_admission_1"
names(studydata)[names(studydata) == "outcome_pat_icu_admission"] <- "Admission_dt_icu_1"
names(studydata)[names(studydata) == "outcome_pat_icu_discharge"] <- "Discharge_dt_icu_1"
names(studydata)[names(studydata) == "outcome_pat_invas_vent"] <- "Invasive_ventilation_1"
names(studydata)[names(studydata) == "outcome_pat_noninvas_vent"] <- "Non_invasive_ventilation_1"
names(studydata)[names(studydata) == "outcome_pat_invas_vent_start"] <- "invasive_start"
names(studydata)[names(studydata) == "outcome_pat_invas_vent_stop"] <- "invasive_stop"
names(studydata)[names(studydata) == "outcome_pat_def_complication.Readmission"] <- "readmissed"

studydata$bacteremia <- ifelse(studydata$`outcome_pat_def_complication#Bacteremia` == 1, "Yes", "No")

studydata$ckd <- ifelse(is.na(studydata$`comorb_pat_def_other#Chronic renal disease`), NA,
                        ifelse(studydata$`comorb_pat_def_other#Chronic renal disease` == 1, "Yes", "No"))
studydata$hypertension <- ifelse(is.na(studydata$`comorb_pat_def_cvd#Hypertension (medicated)`), NA,
                                 ifelse(studydata$`comorb_pat_def_cvd#Hypertension (medicated)` == 1, "Yes", "No"))

## PSI
studydata$`outcome_pat_def_complication#Stroke`

studydata$age_yrs <- as.numeric(studydata$age_yrs)
studydata$PSI_age <- with(studydata, ifelse(gender == "Male", studydata$age_yrs, studydata$age_yrs-10))
studydata$PSI_neoplas <- with(studydata, ifelse(mneoplasm == "Yes", 30, 0)) 
studydata$PSI_liver <- with(studydata, ifelse(`comorb_pat_def_liver#Hepatitis B` == 1 | `comorb_pat_def_liver#Hepatitis C` == 1 | `comorb_pat_def_liver#Cirrhosis` == 1, 20, 0))
studydata$PSI_CHF <- with(studydata, ifelse(`comorb_pat_def_cvd#Congestive heart failure` == 1, 10, 0))
studydata$PSI_CVD <- with(studydata, ifelse(`outcome_pat_def_complication#Stroke` == 1, 10, 0))
studydata$PSI_renal <- with(studydata, ifelse(ckd == "Yes", 10, 0))
studydata$PSI_mental <- with(studydata, ifelse(is.na(base_mental_status), NA, 
                                               ifelse(base_mental_status == 1, 20, 0)))
studydata$PSI_rtr <- with(studydata, ifelse(is.na(rtr), NA,
                                            ifelse(rtr >=30, 20, 0)))                        
studydata$PSI_sys <- with(studydata, ifelse(is.na(sys_bp), NA,
                                            ifelse(sys_bp <=90, 20, 0)))  
studydata$PSI_temp <- with(studydata, ifelse(is.na(Temperature), NA,
                                             ifelse(Temperature <35 | Temperature >=40, 15, 0)))  
studydata$PSI_pulse <- with(studydata, ifelse(is.na(HtR), NA,
                                              ifelse(HtR >=125, 10 , 0)))
studydata$PSI_ABG <- with(studydata, ifelse(is.na(studydata$PH_value_1), NA, 
                                            ifelse(studydata$PH_value_1 <7.35, 30, 0)))
studydata$PSI_urea <- with(studydata, ifelse(!is.na(Blood_Urea_Nitrogen_value_1) & Blood_Urea_Nitrogen_value_1 >= 11 , 20, 0)) ## no NA if missing
studydata$PSI_glucose <- with(studydata, ifelse(!is.na(Glucose_unit_1_1) & Glucose_unit_1_1 >= 14, 10, 0))## no NA if missing
studydata$PSI_Na <- with(studydata, ifelse(is.na(studydata$Sodium_1_1), NA,
                                           ifelse(studydata$Sodium_1_1 <130, 20, 0)))
studydata$PSI_O2 <- with(studydata, ifelse(is.na(PaO2_1) & is.na(oxygen_saturation), NA,
                                           ifelse(!is.na(PaO2_1) & PaO2_1 < 7.99934211, 10,
                                                  ifelse(!is.na(PaO2_1) & oxygen_saturation < 90, 10, 0))))
studydata$PSI_total <- with(studydata, PSI_age + PSI_neoplas + PSI_liver + PSI_CHF + PSI_CVD + PSI_renal + PSI_mental + PSI_rtr + PSI_sys
                            + PSI_temp + PSI_pulse + PSI_ABG + PSI_urea + PSI_Na + PSI_glucose + PSI_O2)
studydata$PSI_calc <- with(studydata, ifelse(PSI_total <=50, 1, 
                                             ifelse(PSI_total %in% 51:70, 2,
                                                    ifelse(PSI_total %in% 71:90, 3,
                                                           ifelse(PSI_total %in% 91:130, 4,
                                                                  ifelse(PSI_total >130, 5, NA))))))

## vanaf EB3.0 niet meer ingevuld maar berekend
studydata$PSI_new <- ifelse(studydata$admission_dt <"2020-10-04", studydata$PSI_score, studydata$PSI_calc)

## antibiotica antibiotic_seven_days
studydata$curmed_pat_antibiotic_start_1 <- as.Date(studydata$curmed_pat_antibiotic_start_1, format = "%d-%m-%Y")
studydata$curmed_pat_antibiotic_start_2 <- as.Date(studydata$curmed_pat_antibiotic_start_2, format = "%d-%m-%Y")
studydata$curmed_pat_antibiotic_start_3 <- as.Date(studydata$curmed_pat_antibiotic_start_3, format = "%d-%m-%Y")
studydata$curmed_pat_antibiotic_start_4 <- as.Date(studydata$curmed_pat_antibiotic_start_4, format = "%d-%m-%Y")
studydata$curmed_pat_antibiotic_start_5 <- as.Date(studydata$curmed_pat_antibiotic_start_5, format = "%d-%m-%Y")
studydata$curmed_pat_antibiotic_start_6 <- as.Date(studydata$curmed_pat_antibiotic_start_6, format = "%d-%m-%Y")

studydata$anti_1 <- as.numeric(difftime(studydata$curmed_pat_antibiotic_start_1, studydata$admission_dt1, units = "days"))
studydata$anti_2 <- as.numeric(difftime(studydata$curmed_pat_antibiotic_start_2, studydata$admission_dt1, units = "days"))
studydata$anti_3 <- as.numeric(difftime(studydata$curmed_pat_antibiotic_start_3, studydata$admission_dt1, units = "days"))
studydata$anti_4 <- as.numeric(difftime(studydata$curmed_pat_antibiotic_start_4, studydata$admission_dt1, units = "days"))
studydata$anti_5 <- as.numeric(difftime(studydata$curmed_pat_antibiotic_start_5, studydata$admission_dt1, units = "days"))
studydata$anti_6 <- as.numeric(difftime(studydata$curmed_pat_antibiotic_start_6, studydata$admission_dt1, units = "days"))

studydata$anti_after_1 <- ifelse(studydata$anti_1 >= 0 & studydata$anti_1 <=7, 1, 0)
studydata$anti_after_2<- ifelse(studydata$anti_2 >= 0 & studydata$anti_2 <=7, 1, 0)
studydata$anti_after_3 <- ifelse(studydata$anti_3 >= 0 & studydata$anti_3 <=7, 1, 0)
studydata$anti_after_4 <- ifelse(studydata$anti_4 >= 0 & studydata$anti_4 <=7, 1, 0)
studydata$anti_after_5 <- ifelse(studydata$anti_5 >= 0 & studydata$anti_5 <=7, 1, 0)
studydata$anti_after_6 <- ifelse(studydata$anti_6 >= 0 & studydata$anti_6 <=7, 1, 0)

antibiotic_7_days <- filter(studydata, studydata$anti_after_1 == 1 | 
                              studydata$anti_after_2 == 1 |
                              studydata$anti_after_3 == 1 |
                              studydata$anti_after_4 == 1 |
                              studydata$anti_after_5 == 1 |
                              studydata$anti_after_6 == 1)
studydata$antibiotic_seven_days <- ifelse(studydata$EB_id %in% antibiotic_7_days$EB_id, "Yes", "No")

## other comorb => all 0
studydata$`comorb_pat_def_immunosup#HIV`
studydata$Dementia <- "No"
studydata$aids_hiv <- ifelse(is.na(studydata$`comorb_pat_def_immunosup#HIV`), NA,
                             ifelse(studydata$`comorb_pat_def_immunosup#HIV`== 1, "Yes", "No"))

## predict variables we want to use
studydata$obesity <- ifelse(studydata$BMI>30, "1", 
                            ifelse(is.na(studydata$BMI), NA, "0"))

## 
studydata$diabetes_comp <- ifelse(studydata$`comorb_pat_organ_fail_dm1#Neuropathy` == 1 | studydata$`comorb_pat_organ_fail_dm1#Nefropathy` == 1 
                                  | studydata$`comorb_pat_organ_fail_dm2#Neuropathy` == 1 | studydata$`comorb_pat_organ_fail_dm2#Nefropathy` == 1
                                  | studydata$`comorb_pat_organ_fail_dm1#Retinopathy` == 1 |studydata$`comorb_pat_organ_fail_dm2#Retinopathy` == 1, 1, 0)

## charlson ## ELDER-BIOME, no peptic ulcer to score, no hemiplegia, rheumatic = connective, mneoplasm has all cancers
## no AIDS only HIV
studydata$metastase <- ifelse(is.na(studydata$metastase), "No", studydata$metastase)
studydata$age <- as.numeric(studydata$age_yrs)
studydata$age_yrs <- round(studydata$age, digits = 0)
studydata$CCI_age <- ifelse(studydata$age <50, 0, ifelse(studydata$age >= 50 & studydata$age < 60, 1,
                                                         ifelse(studydata$age >= 60 & studydata$age < 70, 2,
                                                                ifelse(studydata$age >= 70 & studydata$age < 80, 3,4))))
studydata$cci_mi <- ifelse(studydata$`comorb_pat_def_cvd#Past myocardial infarction`== 1, 1, 0)
studydata$cci_chf <- ifelse(studydata$`comorb_pat_def_cvd#Congestive heart failure` == 1, 1, 0)
studydata$cci_pvd <- ifelse(studydata$`comorb_pat_def_cvd#Peripheral ischaemic/atherosclerotic disease` == "Yes", 1, 0)
studydata$cci_cva <- ifelse(studydata$`outcome_pat_def_complication#Stroke` == 1, 1, 0)
studydata$cci_dementia <- ifelse(studydata$Dementia == "Yes", 1, 0)
studydata$cci_copd <- ifelse(studydata$`comorb_pat_def_cpd#Chronic Obstructive Pulmonary Disease (COPD)` == "Yes", 1, 0)
studydata$cci_ctd <- ifelse(studydata$connective == "Yes", 1, 0)
studydata$cci_pud <- ifelse(studydata$peptic == "Yes", 1, 0)
studydata$cci_liver <- ifelse(studydata$liver_charl == "Mild", 1, ifelse(studydata$liver_charl == "none", 0,3))
studydata$cci_dm <- ifelse(studydata$diabetes == "No" & studydata$diabetes_comp == "No", 0, ifelse(studydata$diabetes_comp == "No", 1, 2))
studydata$cci_hemiplegia <- ifelse(studydata$hemiplegia == "Yes", 2, 0)
studydata$cci_ckd <- ifelse(studydata$ckd == "Yes", 2, 0)
studydata$cci_malig <- ifelse(studydata$metastase == "Yes", 6, 
                              ifelse(studydata$mneoplasm == "Yes", 2, 0))
## studydata$cci_aids <- ifelse(studydata$aids_hiv == "Yes", 6, 0)
studydata$CCI <- studydata$CCI_age + studydata$cci_mi + studydata$cci_chf + studydata$cci_pvd + studydata$cci_cva + studydata$cci_dementia + studydata$cci_copd + studydata$cci_ctd + studydata$cci_pud + studydata$cci_liver + studydata$cci_dm +
  studydata$cci_hemiplegia + studydata$cci_ckd + studydata$cci_malig 
summary(studydata$CCI)

###
studydata$first_admission_dt <- studydata$admission_dt1

## make new variables 
studydata$date_of_death <- as.Date(studydata$date_of_death, format = "%d-%m-%Y")
studydata$Admission_dt_icu_1 <- as.Date(studydata$Admission_dt_icu_1, format = "%d-%m-%Y")
studydata$days_survived_since_ICU <- as.numeric(difftime(studydata$date_of_death, studydata$Admission_dt_icu_1, units = "days"))

studydata$admission_dt1<- as.Date(studydata$admission_dt1, format = "%d-%m-%Y")
studydata$days_survived <- as.numeric(difftime(studydata$date_of_death, studydata$admission_dt1, units = "days"))

## 3wk mortality
studydata$wk3_mortality <- ifelse(is.na(studydata$days_survived), 0, 
                                  ifelse(studydata$days_survived <=21, 1, 0))
## 6wk mortality
studydata$wk6_mortality <- ifelse(is.na(studydata$days_survived), 0, 
                                  ifelse(studydata$days_survived <=42, 1, 0))
## 12wk mortality
studydata$wk12_mortality <- ifelse(is.na(studydata$days_survived), 0, 
                                   ifelse(studydata$days_survived <=84, 1, 0))

studydata$sampling_time <- as.numeric(difftime(studydata$base_date_assessment, studydata$admission_dt1,  units = "days"))
studydata$time_sample_icu <- difftime(studydata$Admission_dt_icu_1, studydata$base_date_assessment, units = "days")

## nog nodig => afname moment biobank
## New variables for both
studydata$Discharge_dt_icu_1 <- as.Date(studydata$Discharge_dt_icu_1, format = "%d-%m-%Y")
studydata$invasive_start <- as.Date(studydata$invasive_start, format = "%d-%m-%Y")
studydata$invasive_stop <- as.Date(studydata$invasive_stop, format = "%d-%m-%Y")
studydata$death_date<- as.Date(studydata$death_date, format = "%d-%m-%Y")
symptom_days <- as.numeric(difftime(studydata$admission_dt1,studydata$onset_dt, units = "days"))
studydata <- add_column(studydata, symptom_days, .after = "admission_dt1")

studydata$discharge_date_local[studydata$MDN == 6271126] <- "2021-05-27"

length_of_stay <- as.numeric(difftime(studydata$discharge_date_local,studydata$admission_dt1, units = "days"))
studydata <- add_column(studydata, length_of_stay, .after = "discharge_date_local")

studydata$oxygen_start <- as.Date(studydata$oxygen_start, format = "%d-%m-%Y")
studydata$oxygen_stop <- as.Date(studydata$oxygen_stop, format = "%d-%m-%Y")
length_of_oxygen <- as.numeric(difftime(studydata$oxygen_stop,studydata$oxygen_start, units = "days"))
studydata <- add_column(studydata, length_of_oxygen, .after = "oxygen_stop")
studydata$day_28_admitted <- ifelse(studydata$length_of_stay >=28, 1, 0)

lenght_of_intubation <- as.numeric(difftime(studydata$invasive_stop,studydata$invasive_start,  units = "days"))
studydata <- add_column(studydata, lenght_of_intubation, .after = "invasive_stop")

## change to numeric 
i <- c("age_yrs", "BMI",	"symptom_days", "CCI",
       "qSOFA_score",	"MEWS_score",	"CURB_score",	"PSI_score",	
       "Haemoglobin_value_1",	"crp_1_1",	"procalcitonine",	"Lactate_2_1",	"LDH",	
       "WBC_2_1",	"Neutrophil_unit_1",	"Lymphocyte_1_1",	
       "Creatinine_value_1",	"Sodium_1_1",	"Potassium_1_1",	"Glucose_unit_1_1",	
       "AST_SGOT_1_1",	"ALT_SGPT_1_1",	"Total_Bilirubin_2_1",	
       "Platelets_value_1",	"d_dimer","pt_spec",	"APT_APTR_1_1",	"PH_value_1",	"PaO2_1",		"PCO2_1",	
       "Blood_Urea_Nitrogen_value_1",	"length_of_stay", "length_of_oxygen")

studydata[ , i] <- apply(studydata[ , i], 2, 
                         function(x) as.numeric(as.character(x)))
remove(i)

studydata$length_of_stay_alive <- ifelse(studydata$hospdeath == 1, NA, studydata$length_of_stay)
studydata$icu_stay <- as.numeric(difftime(studydata$Discharge_dt_icu_1, studydata$Admission_dt_icu_1,  units = "days"))
studydata$days_survived <- difftime(studydata$death_date, studydata$first_admission_dt, unit = "days")
studydata$mortality_d14 <- ifelse(!is.na(studydata$days_survived) & studydata$days_survived <=14, 1, 0)
studydata$mortality_d30 <- ifelse(!is.na(studydata$days_survived) & studydata$days_survived <=30, 1, 0)
studydata$mortality_d90 <- ifelse(!is.na(studydata$days_survived) & studydata$days_survived <=90, 1, 0)


studydata$inclusion_hospital <- ifelse(studydata$patient_volunteer == "Healthy volunteer", "Healthy", studydata$inclusion_hospital)

table(studydata$outcome_pat_cause1)

##
studydata$pathogen_cultured <- ifelse(studydata$outcome_pat_cause1 == 0 | studydata$outcome_pat_cause1 == "No causative agent found",
                                      0, 1)

pneumonia <- filter(studydata, studydata$patient_volunteer == "Pneumonia patient")

## FLOW data
TOM_flow_data <- import("Documents/Luminex/R_code/Original_data/OMIQ_metadata-flow1_metadata name.csv")
pneumonia$flow_Data <- ifelse(pneumonia$EB_id %in% TOM_flow_data$cap_covid_hv, 1, 0)

## spectral flow COVID
ERIK_spectral <- import("Documents/Luminex/R_code/Original_data/spectral data superinfection_from Erik_the COVID-19 patients.xlsx")
pneumonia$spectral <- ifelse(pneumonia$EB_id %in% ERIK_spectral$EB_ID, 1, 0)

##make groups by CAP, COVID-19,hv
#Grouped by pathegon
studydata$covid <- ifelse(studydata$outcome_pat_cause1 == "COVID-19", 1, 
                          ifelse(studydata$outcome_pat_cause2 == "COVID-19", 1, 
                                 ifelse(studydata$outcome_pat_cause3 == "COVID-19", 1, 0)))
studydata$covid <- ifelse(is.na(studydata$covid), 0, studydata$covid)
studydata$COVID_19 <- ifelse(studydata$covid == 1, "Yes", "No")

studydata <- studydata %>%
  mutate(group = case_when(
    EB_id >= 2000 & EB_id < 3000 ~ "HV",
    COVID_19 == "Yes" ~ "COVID",
    COVID_19 == "No" ~ "CAP",
    TRUE ~ NA_character_  # To handle any cases not matching above
  ))
studydata$group

##export form
destination_folder <- "Documents/Luminex/R_code/Original_data/" 
export_file_name <- "studydata.csv" 
write.csv(studydata, file.path(destination_folder, export_file_name), row.names = FALSE)
cat("File exported to:", file.path(destination_folder, export_file_name), "\n")


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
             "hospdeath",	"mortality_d30", "mortality_d90")

catvars <- c("group","flow_Data", "spectral", "inclusion_hospital", "gender",	"ethnic_group#White",	
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
  data        = studydata, 
  strata = "group",
  factorVars  = catvars,
  test        = TRUE)

print(tab1, nonnormal = nonnormal, quote = FALSE, noSpaces = TRUE, smd = T, missing = T)


# Convert the CreateTableOne object to a data frame
tab1_df <- as.data.frame(print(tab1, nonnormal = nonnormal, quote = FALSE, noSpaces = TRUE, smd = TRUE, missing = TRUE))
tab1_df$tablename <- rownames(tab1_df) 
tab1_df <- tab1_df[, c("tablename", setdiff(names(tab1_df), "tablename"))]

# Write the data frame to a CSV file
#export clinical_marker_unique data
destination_folder <- "Documents/Luminex/R_code/Original_data/" 
export_file_name <- "table1.csv" 
write.csv(tab1_df, file.path(destination_folder, export_file_name), row.names = FALSE)
cat("File exported to:", file.path(destination_folder, export_file_name), "\n")




######
one <- studydata$outcome_pat_cause1
two <- studydata$outcome_pat_cause2
three <- studydata$outcome_pat_cause3
all <- c(one, two, three)
print(table(all), quote = T)

##
table(studydata$outcome_pat_cause1)
studydata$strep <- ifelse(studydata$outcome_pat_cause1 == "Streptococcus pneumoniae", 1, 
                          ifelse(studydata$outcome_pat_cause2 == "Streptococcus pneumoniae", 1, 
                                 ifelse(studydata$outcome_pat_cause3 == "Streptococcus pneumoniae", 1, 0)))
studydata$strep <- ifelse(is.na(studydata$strep), 0, studydata$strep)

studydata$covid <- ifelse(studydata$outcome_pat_cause1 == "COVID-19", 1, 
                          ifelse(studydata$outcome_pat_cause2 == "COVID-19", 1, 
                                 ifelse(studydata$outcome_pat_cause3 == "COVID-19", 1, 0)))
studydata$covid <- ifelse(is.na(studydata$covid), 0, studydata$covid)

studydata$aureus <- ifelse(studydata$outcome_pat_cause1 == "Staphylococcus aureus", 1, 
                           ifelse(studydata$outcome_pat_cause2 == "Staphylococcus aureus", 1, 
                                  ifelse(studydata$outcome_pat_cause3 == "Staphylococcus aureus", 1, 0)))
studydata$aureus <- ifelse(is.na(studydata$aureus), 0, studydata$aureus)

studydata$influenza <- ifelse(studydata$outcome_pat_cause1 == "Influenza A virus", 1, 
                              ifelse(studydata$outcome_pat_cause2 == "Influenza A virus", 1, 
                                     ifelse(studydata$outcome_pat_cause3 == "Influenza A virus", 1, 0)))
studydata$influenza <- ifelse(is.na(studydata$influenza), 0, studydata$influenza)

studydata$h_influenz <- ifelse(studydata$outcome_pat_cause1 == "Haemophilus influenzae", 1, 
                               ifelse(studydata$outcome_pat_cause2 == "Haemophilus influenzae", 1, 
                                      ifelse(studydata$outcome_pat_cause3 == "Haemophilus influenzae", 1, 0)))
studydata$h_influenz <- ifelse(is.na(studydata$h_influenz), 0, studydata$h_influenz)

studydata$no_pathogen <- ifelse(studydata$outcome_pat_cause1 == "No causative agent found", 1, 0)

##
studydata$mixed_infections <- ifelse(!is.na(studydata$outcome_pat_cause2), 1, 0)

## 
studydata$other_pathogen <- ifelse(studydata$no_pathogen == 1, 0,
                                   ifelse(studydata$strep == 1, 0, 
                                          ifelse(studydata$covid == 1, 0,
                                                 ifelse(studydata$influenza == 1, 0,
                                                        ifelse(studydata$h_influenz == 1, 0,
                                                               ifelse(studydata$aureus == 1, 0, 1))))))

##
studydata$strep <- ifelse(studydata$mixed_infections == 1, 0, studydata$strep)
studydata$influenza <- ifelse(studydata$mixed_infections == 1, 0, studydata$influenza)
studydata$covid <- ifelse(studydata$mixed_infections == 1, 0, studydata$covid)
studydata$aureus <- ifelse(studydata$mixed_infections == 1, 0, studydata$aureus)
studydata$h_influenz <- ifelse(studydata$mixed_infections == 1, 0, studydata$h_influenz)
studydata$other_pathogen <- ifelse(studydata$mixed_infections == 1, 0, studydata$other_pathogen)

##
other_pathogens <- filter(studydata, studydata$other_pathogen == 1)
other_pathogens <- select(other_pathogens, "EB_id",  "outcome_pat_cause1", "outcome_def_pathogen1")
write.csv(other_pathogens, "other_pathogens.csv")

##
mixed_infection <- filter(studydata, studydata$mixed_infections == 1)
mixed_infection_short <- select(mixed_infection, "EB_id", "outcome_pat_cause1", "outcome_pat_cause2", "outcome_def_pathogen2", "outcome_pat_cause3", 
                                "outcome_def_pathogen3")
write.csv(mixed_infection_short, "mixed_infections.csv")

mixed_infection$outcome_def_pathogen3
mixed_infection$`outcome_pat_microbiology_results_1_Microbiological test 1_Result`
mixed_infection$`outcome_pat_microbiology_results_1_Microbiological test 2_Result`
mixed_infection$`outcome_pat_microbiology_results_1_Microbiological test 3_Result`


##
studydata$pathogen <- ifelse(studydata$mixed_infections == 1, "mixed_infection", 
                             ifelse(studydata$no_pathogen == 1, "No_pathogen", 
                                    ifelse(studydata$strep == 1, "Strep", 
                                           ifelse(studydata$aureus == 1, "aureus", 
                                                  ifelse(studydata$influenza == 1, "influenza", 
                                                         ifelse(studydata$h_influenz == 1, "haemo", 
                                                                ifelse(studydata$covid == 1, "covid", "other")))))))
table(studydata$pathogen)



##



studydata$inclusion_hospital

##
studydata$outcome_pat_cause1 <- ifelse(is.na(studydata$outcome_pat_cause1), 0, studydata$outcome_pat_cause1)
studydata$outcome_pat_cause2 <- ifelse(is.na( studydata$outcome_pat_cause2), 0, studydata$outcome_pat_cause2)
studydata$outcome_pat_cause3 <- ifelse(is.na(studydata$outcome_pat_cause3), 0, studydata$outcome_pat_cause3)

## "Corona virus" , "Stenotrophomonas maltophilia", "Moraxella osloensis",  "Streptococcus salivarius"
studydata$Fungi <- ifelse(studydata$outcome_pat_cause1 ==  "Aspergillus spp" | studydata$outcome_pat_cause1 ==  "Pneumocystis carinii", 1, 0)
studydata$Gram_negative <- ifelse(studydata$outcome_pat_cause1 ==  "Escherichia coli" | studydata$outcome_pat_cause1 == "Haemophilus influenzae" |  studydata$outcome_pat_cause1 == "Pseudomonas aeruginosa" , 1, 0)
studydata$Gram_positive <- ifelse(studydata$outcome_pat_cause1 ==  "Rothia dentocariosa" | studydata$outcome_pat_cause1 ==  "Staphylococcus aureus" | studydata$outcome_pat_cause1 ==  "Streptococcus pneumoniae" 
                                  | studydata$outcome_pat_cause1 ==  "Streptococcus pyogenes", 1, 0)
studydata$Virus <- ifelse(studydata$outcome_pat_cause1 ==  "Coronavirus" | studydata$outcome_pat_cause1 == "human MetaPneumoVirus (hMPV)" | studydata$outcome_pat_cause1 == "Influenza A virus" 
                          |  studydata$outcome_pat_cause1 == "Influenza B virus" | studydata$outcome_pat_cause1 == "Parainfluenza virus 1-4" | studydata$outcome_pat_cause1 == "Respiratory syncytial virus (RSV)"
                          |  studydata$outcome_pat_cause1 == "Rhinovirus", 1, 0)
studydata$Other <- ifelse(studydata$outcome_pat_cause1 ==  "Mycobacterium tuberculosis", 1, 0)
studydata$Unknown <- ifelse(studydata$outcome_pat_cause1 ==  "No causative agent found", 1, 0)

studydata$Fungi <- ifelse(studydata$outcome_pat_cause2 ==  "Aspergillus spp" | studydata$outcome_pat_cause2 ==  "Pneumocystis carinii", 1, studydata$Fungi)
studydata$Gram_negative <- ifelse(studydata$outcome_pat_cause2 ==  "Escherichia coli" | studydata$outcome_pat_cause2 == "Haemophilus influenzae" |  studydata$outcome_pat_cause2 == "Pseudomonas aeruginosa" , 1, 
                                  studydata$Gram_negative)
studydata$Gram_positive <- ifelse(studydata$outcome_pat_cause2 ==  "Rothia dentocariosa" | studydata$outcome_pat_cause2 ==  "Staphylococcus aureus" | studydata$outcome_pat_cause2 ==  "Streptococcus pneumoniae" 
                                  | studydata$outcome_pat_cause2 ==  "Streptococcus pyogenes", 1, studydata$Gram_positive)
studydata$Virus <- ifelse(studydata$outcome_pat_cause2 ==  "Coronavirus" | studydata$outcome_pat_cause2 == "human MetaPneumoVirus (hMPV)" | studydata$outcome_pat_cause2 == "Influenza A virus" 
                          |  studydata$outcome_pat_cause2 == "Influenza B virus" | studydata$outcome_pat_cause2 == "Parainfluenza virus 1-4" | studydata$outcome_pat_cause2 == "Respiratory syncytial virus (RSV)"
                          |  studydata$outcome_pat_cause2 == "Rhinovirus", 1, studydata$Virus)
studydata$Other <- ifelse(studydata$outcome_pat_cause2 ==  "Mycobacterium tuberculosis", 1, studydata$Other)

studydata$Fungi <- ifelse(studydata$outcome_pat_cause3 ==  "Aspergillus spp" | studydata$outcome_pat_cause3 ==  "Pneumocystis carinii", 1, studydata$Fungi)
studydata$Gram_negative <- ifelse(studydata$outcome_pat_cause3 ==  "Escherichia coli" | studydata$outcome_pat_cause3 == "Haemophilus influenzae" |  studydata$outcome_pat_cause3 == "Pseudomonas aeruginosa" , 1, 
                                  studydata$Gram_negative)
studydata$Gram_positive <- ifelse(studydata$outcome_pat_cause3 ==  "Rothia dentocariosa" | studydata$outcome_pat_cause3 ==  "Staphylococcus aureus" | studydata$outcome_pat_cause3 ==  "Streptococcus pneumoniae" 
                                  | studydata$outcome_pat_cause3 ==  "Streptococcus pyogenes", 1, studydata$Gram_positive)
studydata$Virus <- ifelse(studydata$outcome_pat_cause3 ==  "Coronavirus" | studydata$outcome_pat_cause3 == "human MetaPneumoVirus (hMPV)" | studydata$outcome_pat_cause3 == "Influenza A virus" 
                          |  studydata$outcome_pat_cause3 == "Influenza B virus" | studydata$outcome_pat_cause3 == "Parainfluenza virus 1-4" | studydata$outcome_pat_cause3 == "Respiratory syncytial virus (RSV)"
                          |  studydata$outcome_pat_cause3 == "Rhinovirus", 1, studydata$Virus)
studydata$Other <- ifelse(studydata$outcome_pat_cause3 ==  "Mycobacterium tuberculosis", 1, studydata$Other)

##
studydata$outcome_def_pathogen2 <- ifelse(is.na(studydata$outcome_def_pathogen2), 0, studydata$outcome_def_pathogen2)
studydata$outcome_def_pathogen1 <- ifelse(is.na( studydata$outcome_def_pathogen1), 0, studydata$outcome_def_pathogen1)

studydata$Gram_negative <- ifelse(studydata$outcome_def_pathogen2 == "Moraxella osloensis" | studydata$outcome_def_pathogen2 == "Stenotrophomonas maltophilia", 1, studydata$Gram_negative)
studydata$Gram_positive <- ifelse(studydata$outcome_def_pathogen1 == "Streptococcus salivarius", 1, studydata$Gram_positive)
studydata$Virus <- ifelse(studydata$outcome_def_pathogen1 == "Corona virus", 1, studydata$Virus)

studydata$outcome_pat_cause1



#Question:
#need to check the mortality with code from EB
#imputaed MEWS?
#if the ID not retested for COVID-19 is right?