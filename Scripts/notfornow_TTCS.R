#####                #####
###   CALCULATE TTCS   ###
#####                #####


#getwd()

#pneumonia <- filter(sanquin, sanquin$patient_volunteer == "Pneumonia patient")
pneumonia <- read.csv("Original_data/Luminex_CAP_COVID_HV_log.csv")

df_tcs <- import("Original_data/ELDER-BIOME_excel_export_20230418012733.xlsx")
## df_tcs <- df_tcs[df_tcs$`Participant Id` %in% pneumonia$EB_id, ] ### REPLACE WITH YOUR DF
n_distinct(df_tcs$`Participant Id`) ## 839 => correct

openxlsx::write.xlsx(df_tcs, "TTCS_raw.xlsx")

##
pneumonia$Record_ID <- pneumonia$ID
df_tcs$Record_ID <- df_tcs$`Participant Id`

df_tcs <- left_join(df_tcs, pneumonia[,c('Record_ID', 'admission_dt1', 'discharge_date_local', 'day28_pat_status', 'date_of_death')], by = 'Record_ID')
df_tcs$clin_stab_date <- as.Date(df_tcs$clin_stab_date, format = "%d-%m-%Y")
df_tcs$los <- as.numeric(difftime(df_tcs$discharge_date_local, df_tcs$admission_dt1, unit = 'days'))


## original HALM temp <= 37.2, htr <100, rr <=24, sysbp >=90, saturation >=90% of sAo2 >=60 on room air
## modified HALM temp <= 37.8
df_tcs$cs_temp <- with(df_tcs, ifelse(is.na(clin_stab_temp), 1,
                                      ifelse(clin_stab_temp <= 37.2, 1, 0)))
df_tcs$cs_hr <- with(df_tcs, ifelse(is.na(clin_stab_hr), 1,
                                    ifelse(clin_stab_hr <= 100, 1, 0)))
df_tcs$cs_rr <- with(df_tcs, ifelse(is.na(clin_stab_resp), 1,
                                    ifelse(clin_stab_resp <= 24, 1, 0)))
df_tcs$cs_bp <- with(df_tcs, ifelse(is.na(clin_stab_sbp), 1,
                                    ifelse(clin_stab_sbp >= 90, 1, 0)))
df_tcs$cs_ox <- with(df_tcs, ifelse(clin_stab_add_oxygen == "No" & clin_stab_oxygen >=90, 1, 0))

##
df_tcs$stable <- with(df_tcs, ifelse(cs_temp == 1 & cs_hr == 1 & cs_rr == 1 & cs_bp == 1 & cs_ox == 1, 1, 0))
df_tcs$stable[df_tcs$clin_stab_hosp == 0] <- 0 # you can't be stable upon admission

df_stable <- df_tcs %>% arrange(Record_ID, desc(stable), clin_stab_hosp) %>% distinct(Record_ID, .keep_all = T)
df_stable$clin_stab_hosp_old <- ifelse(df_stable$clin_stab_hosp == 0 & df_stable$stable == 0 & df_stable$day28_pat_status == 5, 28,
                                   ifelse(df_stable$clin_stab_hosp == 0 & df_stable$stable == 0, df_stable$los, df_stable$clin_stab_hosp))

## original HALM temp <= 37.2, htr <100, rr <=24, sysbp >=90, saturation >=90% of sAo2 >=60 on room air
## modified HALM temp <= 37.8
df_tcs$cs_temp <- with(df_tcs, ifelse(is.na(clin_stab_temp), 1,
                                      ifelse(clin_stab_temp <= 37.8, 1, 0)))
df_tcs$cs_hr <- with(df_tcs, ifelse(is.na(clin_stab_hr), 1,
                                    ifelse(clin_stab_hr <= 100, 1, 0)))
df_tcs$cs_rr <- with(df_tcs, ifelse(is.na(clin_stab_resp), 1,
                                    ifelse(clin_stab_resp <= 24, 1, 0)))
df_tcs$cs_bp <- with(df_tcs, ifelse(is.na(clin_stab_sbp), 1,
                                    ifelse(clin_stab_sbp >= 90, 1, 0)))

df_tcs$cs_ox <- with(df_tcs, ifelse(!is.na(clin_stab_add_oxygen) & clin_stab_add_oxygen == "No" & clin_stab_oxygen >=90, 1, 0))
table(df_tcs$clin_stab_add_oxygen, df_tcs$cs_ox)

##
df_tcs$stable <- with(df_tcs, ifelse(cs_temp == 1 & cs_hr == 1 & cs_rr == 1 & cs_bp == 1 & cs_ox == 1, 1, 0))
df_tcs$stable[df_tcs$clin_stab_hosp == 0] <- 0 # you can't be stable upon admission

df_stable_long <- df_tcs %>% arrange(Record_ID, desc(stable), clin_stab_hosp) %>% distinct(Record_ID, .keep_all = T)
df_stable_long$clin_stab_hosp_new <- ifelse(df_stable_long$clin_stab_hosp == 0 & df_stable_long$stable == 0 & df_stable_long$day28_pat_status == 5, 28,
                                         ifelse(df_stable_long$clin_stab_hosp == 0 & df_stable_long$stable == 0, df_stable_long$los, df_stable_long$clin_stab_hosp))

## combine
df_stable_long <- select(df_stable_long, "Record_ID", "clin_stab_hosp_new")
df_stable <- merge(df_stable, df_stable_long, by = "Record_ID")
df_stable <- select(df_stable, "Record_ID", "clin_stab_hosp_old", "clin_stab_hosp_new")
openxlsx::write.xlsx(df_stable, "TTCS_cleaned.xlsx")

##
sanquin <- merge(sanquin, df_stable, by.x = "EB_id", by.y = "Record_ID", all.x = T)

df_stable$Record_ID
