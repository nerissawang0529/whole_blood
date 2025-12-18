#rm(list = ls())

## install.packages("pacman")
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
               tableone, kmed, survival, survminer, skimr, lattice, reshape2, grDevices, ggrepel)


#Ready for luminex data
Luminex_CAP_COVID_HV <- read.csv("Original_data/Luminex_CAP_COVID_HV_raw.csv")


#Ready for studydata (copy the code for studydata for HV and patients)
studydata <- import("Original_data/ELDER-BIOME_excel_export_20230418012733.xlsx")
#this is the raw data for study
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

#remove the ones not retesed (COVID-19)
ids_to_remove <- c(3024, 1194, 1195, 1210, 1227, 1228, 3095, 3098, 3099, 3104, 
                   3115, 3128, 3165, 3166, 3167, 3181, 3186, 3201, 3207, 3208, 
                   3218, 3219, 3227, 3245, 3286, 1149, 1160, 1186, 3174, 3291, 3248)

studydata <- studydata %>%
  filter(!`Participant Id` %in% ids_to_remove)

#exclude which don't have M value (too high and not retested)
ids_to_remove_2 <- c(1069,1004)
studydata <- studydata %>%
  filter(!`Participant Id` %in% ids_to_remove_2)

colnames(studydata)[colnames(studydata) == "Participant Id"] <- "Record Id"
names(studydata)[names(studydata) == "Record Id"] <- "EB_id"


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

notCOVID_notCAP <- studydata[studydata$group == "CAP" & studydata$outcome_pat_posthoc == "It was not a pneumonia", ]
notCOVID_notCAP <- notCOVID_notCAP[,c("EB_id","outcome_pat_posthoc","outcome_pat_posthoc_other","outcome_pat_other_diagnosis")]

#studydata data for CAP
studydata <- studydata %>%
  filter(!EB_id %in% notCOVID_notCAP$EB_id)


#we want to figure out if the biopeak is because different time batch (different person who collected at different time)
studydata$date_created <- as.POSIXct(strptime(studydata$`Participant Creation Date`, format="%d-%m-%Y %H:%M:%S"))

table(studydata$group)

datab <- merge(Luminex_CAP_COVID_HV, studydata, 
              by.x = "ID", by.y = "EB_id")

table(datab$group)

datab$date_created2 <- as.Date(datab$date_created)

#draw the time how many patients were collected at different time
ggplot(data = datab, aes(x = date_created2)) +
  geom_histogram(binwidth = 30) +  # Counts occurrences per date
  facet_wrap(~ group, nrow=3) +
  scale_x_date(date_labels = "%Y-%m-%d", date_breaks = "2 weeks") 

#double check the biopeak for patients
hist(log10(datab$IL_8_1[datab$stimulation=="M"]))
hist(log10(datab$IL_8_1[datab$stimulation=="M" & datab$group=="CAP" ]))
hist(log10(datab$IL_8_1[datab$stimulation=="M" & datab$group=="HV" ]))


#we only interested in medium
databM <- datab[datab$stimulation=="M",]
#only check the day1 and day2
#databM <- databM[databM$Day %in% c("Day_1", "Day_2"), ] #we are interested in all data, so not run this line

#we choosed IL_6 to try to make the patients into groups (we want to match the high IL_6 into the time).
databM$highIL6 <- log10(databM$IL_6) > 3
table(databM$highIL6 )

#group the patients in HV AND CAP
HV_highIL6 <- unique( databM$ID[log10(databM$IL_6) > 3 & databM$group=="HV"] )
CAP_highIL6 <- unique( databM$ID[log10(databM$IL_6) > 3 & databM$group=="CAP"] )


#Check the High IL6, if they are in one batch time?
ggplot(data=databM, aes(x=as.factor(highIL6), y=date_created2))+
  geom_boxplot()+
  facet_wrap(~ group, ncol=3)+
  stat_compare_means(method = "t.test", label = "p.format")
#the answer is no.

#make the table one for high IL6 AND low IL6 groups
##HV
studydata_HV <- read.csv("Original_data/studydata_HV.csv")

studydata_HV$highIL6 <- studydata_HV$EB_id %in% HV_highIL6
table(studydata_HV$highIL6)

## table
allvars <- c("inclusion_hospital", "age_yrs",	"gender",	"ethnic_group#White",	
             "BMI",	"symptom_days",
             "CCI",	"hypertension",	"cpd",	"COPD",	"diabetes",	
             "ccd","ckd",	"mneoplasm",	"immune_sup",	"cnd",	
             "qSOFA_score",	"MEWS_score",	"CURB_score",	"PSI_new",
             "crp_1_1",	 "WBC_2_1",	"Neutrophil_unit_1",	"Lymphocyte_1_1",	
             "Creatinine_value_1",
             "Platelets_value_1",	
             "Blood_Urea_Nitrogen_value_1",	"antibiotic_seven_days", 
             "length_of_stay", 
             "hospdeath",	"mortality_d30", "mortality_d90","highIL6")

catvars <- c("inclusion_hospital", "gender",	"ethnic_group#White",	
             "hypertension",	"cpd",	"COPD",	"diabetes",	
             "ccd","ckd",	"mneoplasm",	"immune_sup",	"cnd",
             "antibiotic_seven_days", 
             "hospdeath",	"mortality_d30", "mortality_d90","highIL6")

nonnormal <- c("age_yrs", "BMI", "symptom_days",  "CCI", "qSOFA_score",	"MEWS_score",	"CURB_score",	"PSI_new",
               "crp_1_1",	 "WBC_2_1",	"Neutrophil_unit_1",	"Lymphocyte_1_1",	
               "Creatinine_value_1",
               "Platelets_value_1",	
               "Blood_Urea_Nitrogen_value_1",	"length_of_oxygen","length_of_stay" ,"icu_stay")

tab2     <- CreateTableOne(
  vars        = allvars, 
  data        = studydata_HV, 
  strata = "highIL6",
  factorVars  = catvars,
  test        = TRUE)

print(tab2, nonnormal = nonnormal, quote = TRUE, noSpaces = TRUE, smd = T, missing = T)

#CAP####
studydata_patients <- read.csv("Original_data/studydata_patients.csv")
studydata_patients$highIL6 <- studydata_patients$EB_id %in% CAP_highIL6
## table
allvars <- c("inclusion_hospital", "age_yrs",	"gender",	"ethnic_group#White",	
             "BMI",	"symptom_days",
             "CCI",	"hypertension",	"cpd",	"COPD",	"diabetes",	
             "ccd","ckd",	"mneoplasm",	"immune_sup",	"cnd",	
             "qSOFA_score",	"MEWS_score",	"CURB_score",	"PSI_new",
             "crp_1_1",	 "WBC_2_1",	"Neutrophil_unit_1",	"Lymphocyte_1_1",	
             "Creatinine_value_1",
             "Platelets_value_1",	
             "Blood_Urea_Nitrogen_value_1",	"antibiotic_seven_days", 
             "length_of_stay", 
             "hospdeath",	"mortality_d30", "mortality_d90","highIL6")

catvars <- c("inclusion_hospital", "gender",	"ethnic_group#White",	
             "hypertension",	"cpd",	"COPD",	"diabetes",	
             "ccd","ckd",	"mneoplasm",	"immune_sup",	"cnd",
             "antibiotic_seven_days", 
             "hospdeath",	"mortality_d30", "mortality_d90","highIL6")

nonnormal <- c("age_yrs", "BMI", "symptom_days",  "CCI", "qSOFA_score",	"MEWS_score",	"CURB_score",	"PSI_new",
               "crp_1_1",	 "WBC_2_1",	"Neutrophil_unit_1",	"Lymphocyte_1_1",	
               "Creatinine_value_1",
               "Platelets_value_1",	
               "Blood_Urea_Nitrogen_value_1",	"length_of_oxygen","length_of_stay" ,"icu_stay")

tab2     <- CreateTableOne(
  vars        = allvars, 
  data        = studydata_patients, 
  strata = "highIL6",
  factorVars  = catvars,
  test        = TRUE)

print(tab2, nonnormal = nonnormal, quote = TRUE, noSpaces = TRUE, smd = T, missing = T)


#if the high medium has higher stimulation?
merged_data <- merge(Luminex_CAP_COVID_HV, studydata_patients, 
                     by.x = "ID", by.y = "EB_id")
merged_data_M <-  merged_data[merged_data$stimulation == "M",]
plot(log10(merged_data_M$crp_1_1) , log10(merged_data_M$IL_8_1 ))


hist(log10(databM$IL_6))


#mimic the figures from Xanthe Brands
ggplot(data=datab, aes(x= as.factor(group), y=log10(IL_6) ))+
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.1, binwidth = 0.05)

ggplot(data=datab, aes(x=log10(IL_6) ))+
  geom_histogram(aes(y=..density..), alpha=0.8, binwidth = 0.2)+
  facet_wrap(~group)+
  theme_classic()


ggplot(data=datab, aes(x=log10(IL_6) ))+
  geom_histogram(aes(y=..density..), alpha=0.8, binwidth = 0.2)+
  facet_wrap(~group * stimulation)+
  theme_classic()


datab$BasXan <- datab$date_created2 < "2018-07-30"


databBX <- datab[datab$BasXan & !(log10(datab$IL_6) < 0.3 & 
                                    datab$group=="HV" & 
                                    datab$stimulation=="Kpneu"),]

ggplot(data=databBX, aes(x=log10(IL_6) ))+
  geom_histogram(aes(y=..density..), alpha=0.8, binwidth = 0.05)+
  facet_wrap(~group * stimulation , scales = "free")+
  theme_classic()



databBX <- datab[datab$BasXan & !(datab$ID==2005),]


HVBX <- databBX[databBX$group=="HV",]

plot( log10(HVBX$IL_6[HVBX$stimulation=="LPS"]),
      log10(HVBX$IL_6[HVBX$stimulation=="Kpneu"])      
)

library(reshape2)
mdataBX <- HVBX[,1:12]
mdataBX$Day <- NULL

ldataBX <- reshape2::melt(mdataBX, id=c("ID","stimulation") )
ldataBX$level <- log10(ldataBX$value)

# Vertical box plot by group
boxplot(level ~ stimulation, data = ldataBX[ldataBX$variable=="IL_6" & 
                                              ldataBX$stimulation %in% c("Kpneu","LPS"),], 
        col = "white")

# Points
stripchart(level ~ stimulation,
           data = ldataBX[ldataBX$variable=="IL_6" & ldataBX$stimulation %in% c("Kpneu","LPS") ,],
           method = "jitter",
           pch = 19,
           col = 2:4,
           vertical = TRUE,
           add = TRUE)


## highest 9 LPS

LPSBX <- ldataBX[ldataBX$variable=="IL_6" & 
                   ldataBX$stimulation %in% c("LPS"),]

top9 <- head(LPSBX[order(-LPSBX$value ),], 9)

HVBX$top9 <- as.character(HVBX$ID) %in% as.character(top9$ID)

table(HVBX$top9)

plot( log10(HVBX$IL_6[HVBX$stimulation=="LPS"]),
      log10(HVBX$IL_6[HVBX$stimulation=="Kpneu"]) ,
      col= as.factor(HVBX$top9[HVBX$stimulation=="LPS"]),
      pch=16
)


#IL_6 VS time #to check if is because somthing wrong with time
library(ggplot2)
library(dplyr)

#ONLY HV
library(ggplot2)
library(dplyr)

# 
databM_HV <- databM %>% filter(group == "HV")

# Step 1：IL-6
ggplot(data = databM_HV, aes(x = date_created2, y = log10(IL_6))) +
  geom_point(color = "forestgreen", alpha = 0.7, size = 1.8) +
  scale_x_date(date_breaks = "3 month", date_labels = "%Y-%m") +
  labs(x = "Sample Collection Date (by 3 Months)", y = "log10(IL-6)", title = "IL-6 levels in non-infection group under medium over time") +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
library(tidyr)

# all marker
markers <- c("IL_8_1", "IL_6", "IL_10", "IL_1RA", "IL_1beta", "TNF", "CCL2", "CCL3", "CCL4")

# 
datab_long_HV <- databM_HV %>%
  pivot_longer(cols = all_of(markers), names_to = "marker", values_to = "value") %>%
  mutate(log_value = log10(value))

# Step 2
ggplot(datab_long_HV, aes(x = date_created2, y = log_value)) +
  geom_point(color = "forestgreen", alpha = 0.7, size = 1.5) +
  facet_wrap(~ marker, scales = "free_y", ncol = 3) +
  scale_x_date(date_breaks = "1 month", date_labels = "%Y-%m") +
  labs(x = "Sample Collection Date (by Month)", y = "log10(Marker Value)", title = "Cytokine markers in HV group over time") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




#HV AND CAP
databM_filtered <- databM %>% filter(group != "COVID")
#only IL_6
ggplot(data = databM_filtered, aes(x = date_created2, y = log10(IL_6), color = group)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_x_date(date_breaks = "1 month", date_labels = "%Y-%m") +
  labs(x = "Sample Collection Date (by Month)", y = "log10(IL-6)", color = "Group") +
  theme_classic(base_size = 12) +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#all markers
library(dplyr)
library(ggplot2)
library(tidyr)

databM_filtered <- databM %>%
  filter(group != "COVID")

markers <- c("IL_8_1", "IL_6", "IL_10", "IL_1RA", "IL_1beta", "TNF", "CCL2", "CCL3", "CCL4")

datab_long <- databM_filtered %>%
  pivot_longer(cols = all_of(markers), names_to = "marker", values_to = "value") %>%
  mutate(log_value = log10(value))

ggplot(datab_long, aes(x = date_created2, y = log_value, color = group)) +
  geom_point(alpha = 0.6, size = 1.5) +
  facet_wrap(~ marker, scales = "free_y", ncol = 3) +
  scale_x_date(date_breaks = "1 month", date_labels = "%Y-%m") +
  labs(x = "Sample Collection Date (by Month)", y = "log10(Marker Value)", color = "Group") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


highIL6 <- c(HV_highIL6, CAP_highIL6)

datab_long$HiIL6 <- ifelse(datab_long$ID %in% highIL6, 
                           "HighIL6_M", "LowIL6_M")




plot5 <- ggplot(datab_long, aes(x = date_created2, y = log_value, color = HiIL6)) +
  geom_point(alpha = 0.6, size = 1.5) +
  facet_grid(marker ~ group, scales = "free_y") +
  scale_x_date(date_breaks = "3 month", date_labels = "%Y-%m") +
  labs(x = "Sample Collection Date (by Month)", y = "log10(Marker Value)", color = "Group") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot5
ggsave(file="/Users/huiwang/Documents/Luminex/luminex_R/Figure/months.svg", plot=plot5, width=15, height=12)




# grouped by IL_6 to check if they are the same group of hv for other markers
#
library(ggplot2)
library(dplyr)

# 
hv_il6 <- datab %>%
  filter(stimulation == "M", group == "HV", !is.na(IL_6)) %>%
  mutate(log_IL6 = log10(IL_6))

# 
cutoff_q3 <- quantile(hv_il6$log_IL6, 0.75, na.rm = TRUE)

# 
ggplot(hv_il6, aes(x = log_IL6)) +
  geom_density(fill = "#A9C1D9", color = "#2A4765", alpha = 0.85, linewidth = 1.2) +
  #geom_vline(xintercept = cutoff_q3, color = "#1D3557", linetype = "dashed", linewidth = 1) +
  labs(
    x = expression(log[10]~"(IL-6)"),
    y = "Density",
    title = "IL-6 Distribution (Medium Stimulation, Non-Infection Group)",
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

# 设置 markers
markers <- c("IL_8_1", "IL_10", "IL_1RA", "IL_1beta", "TNF", "CCL2", "CCL3", "CCL4")

# 创建 IL-6 分组变量（基于 HV 的 Q3 cutoff）
cutoff_q3 <- quantile(log10(datab$IL_6[datab$stimulation == "M" & datab$group == "HV"]), 0.75, na.rm = TRUE)
datab$highIL6 <- log10(datab$IL_6) > cutoff_q3

# 整理数据（HV + medium 刺激）
hv_filtered <- datab %>%
  filter(group == "HV", stimulation == "M") %>%
  dplyr::select(ID, all_of(markers), highIL6) %>%
  pivot_longer(cols = all_of(markers), names_to = "marker", values_to = "value") %>%
  mutate(log_value = log10(value),
         highIL6 = factor(highIL6, levels = c(FALSE, TRUE), labels = c("Low IL-6", "High IL-6")))

# 绘图
ggplot(hv_filtered, aes(x = highIL6, y = log_value, fill = highIL6)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8, color = "black") +
  geom_jitter(width = 0.15, size = 1.3, alpha = 0.5, color = "black") +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 3.5) +
  facet_wrap(~ marker, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c("Low IL-6" = "green3", "High IL-6" = "red")) +
  labs(
    x = "",
    y = "log10(Marker Value)",
    title = "Marker Levels by IL-6 Group (Full Erik & Nerissa Dataset)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(face = "bold")
  )

#if the high part higher after stimulation?
library(ggplot2)
library(dplyr)
library(tidyr)

plot_df <- datab_HV %>%
  dplyr::select(ID, stimulation, log_IL8, highIL6) %>%
  distinct() %>%
  pivot_wider(names_from = stimulation, values_from = log_IL8) %>%
  filter(!is.na(M), !is.na(LPS))  

plot_df_long <- plot_df %>%
  pivot_longer(cols = c("M", "LPS"), names_to = "stimulation", values_to = "log_IL8")

ggplot(plot_df_long, aes(x = stimulation, y = log_IL8, group = ID)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_line(aes(color = stimulation), size = 1.2, alpha = 0.6) +
  facet_wrap(~ highIL6, labeller = as_labeller(c(`TRUE` = "High IL-6", `FALSE` = "Low IL-6"))) +
  scale_color_manual(values = c("M" = "#999999", "LPS" = "#1f77b4")) +
  labs(
    title = "IL-8 Change (Medium vs LPS) by IL-6 Group",
    x = NULL,
    y = "log10(IL-8)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

plot_df %>%
  mutate(delta = LPS - M,
         IL6_group = ifelse(highIL6, "High IL-6", "Low IL-6")) %>%
  ggplot(aes(x = IL6_group, y = delta, fill = IL6_group)) +
  geom_boxplot(alpha = 0.5, width = 0.4, color = "black") +
  geom_jitter(width = 0.15, alpha = 0.8, size = 1.5) +
  scale_fill_manual(values = c("High IL-6" = "#003f5c", "Low IL-6" = "#7a5195")) +
  labs(
    x = NULL, y = "IL-8 Response (LPS - Medium)", 
    title = "Delta IL-8 Response by IL-6 Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5)) +
  stat_compare_means(method = "t.test", label = "p.format")


#for three different stimulation
library(ggplot2)
library(dplyr)
library(tidyr)

datab_HV <- datab %>% filter(group == "HV", stimulation %in% c("M", "LPS", "Kpneu", "Spneu"))

cutoff_q3 <- quantile(log10(datab$IL_6[datab$stimulation == "M" & datab$group == "HV"]), 0.75, na.rm = TRUE)

datab_HV <- datab_HV %>%
  group_by(ID) %>%
  mutate(highIL6 = any(log10(IL_6[stimulation == "M"]) > cutoff_q3)) %>%
  ungroup()

plot_df_long <- datab_HV %>%
  dplyr::select(ID, stimulation, IL_8_1, highIL6) %>%
  mutate(log_IL8 = log10(IL_8_1)) %>%
  filter(!is.na(log_IL8))

plot_df_long$stimulation <- factor(plot_df_long$stimulation, levels = c("M", "LPS", "Kpneu", "Spneu"))

ggplot(plot_df_long, aes(x = stimulation, y = log_IL8, group = ID)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_line(aes(color = stimulation), size = 1.2, alpha = 0.6) +
  facet_wrap(~ highIL6, labeller = as_labeller(c(`TRUE` = "High IL-6", `FALSE` = "Low IL-6"))) +
  scale_color_manual(values = c("M" = "#999999", "LPS" = "#1f77b4", "Kpneu" = "#e377c2", "Spneu" = "#2ca02c")) +
  labs(
    title = "IL-8 Response Across Stimulations by IL-6 Baseline",
    x = NULL,
    y = "log10(IL-8)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

markers_to_plot <- c("IL_8_1", "IL_10", "TNF", "CCL2", "CCL3", "CCL4", "IL_1RA", "IL_1beta")
library(dplyr)
library(tidyr)

# 仅选 HV 组 + 目标刺激
datab_HV <- datab %>% 
  filter(group == "HV", stimulation %in% c("M", "LPS", "Kpneu", "Spneu"))

# 计算 IL-6 Q3 cutoff，用于 high/low IL-6 baseline 分组
cutoff_q3 <- quantile(log10(datab$IL_6[datab$stimulation == "M" & datab$group == "HV"]), 0.75, na.rm = TRUE)

# 标注 highIL6 分组
datab_HV <- datab_HV %>%
  group_by(ID) %>%
  dplyr::mutate(highIL6 = any(log10(IL_6[stimulation == "M"]) > cutoff_q3)) %>%
  ungroup()

# pivot 为 long format，包含 log10 值
long_markers <- datab_HV %>%
  dplyr::select(ID, stimulation, highIL6, all_of(markers_to_plot)) %>%
  pivot_longer(cols = all_of(markers_to_plot), names_to = "marker", values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(log_value = log10(value)) %>%
  filter(!is.infinite(log_value))  # 去除 log10(0) 的 -Inf

long_markers$stimulation <- factor(long_markers$stimulation, levels = c("M", "LPS", "Kpneu", "Spneu"))
ggplot(long_markers, aes(x = stimulation, y = log_value, group = ID)) +
  geom_point(size = 1.5, alpha = 0.7) +
  geom_line(aes(color = stimulation), size = 1.1, alpha = 0.5) +
  facet_wrap(marker ~ highIL6, ncol = 6, labeller = labeller(
    highIL6 = c(`TRUE` = "High IL-6", `FALSE` = "Low IL-6")
  )) +
  scale_color_manual(values = c("M" = "#999999", "LPS" = "#1f77b4", "Kpneu" = "#e377c2", "Spneu" = "#2ca02c")) +
  labs(
    title = "Marker Responses Across Stimulations (Non-Infection Patients)",
    x = NULL,
    y = "log10(Marker Value)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )
ggplot(long_markers, aes(x = stimulation, y = log_value, group = ID)) +
  geom_point(aes(color = highIL6), size = 1.5, alpha = 0.8) +
  geom_line(aes(color = highIL6), size = 0.2, alpha = 0.5) +
  facet_wrap(marker ~ highIL6, ncol = 8, labeller = labeller(
    highIL6 = c(`TRUE` = "High IL-6", `FALSE` = "Low IL-6")
  )) +
  scale_color_manual(
    values = c(`FALSE` = "green3", `TRUE` = "red"),
    labels = c("Low IL-6", "High IL-6")
  ) +
  labs(
    title = "Marker Responses Across Stimulations (Non-Infection Patients)",
    x = NULL,
    y = "log10(Marker Value)",
    color = "IL-6 Group"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )


