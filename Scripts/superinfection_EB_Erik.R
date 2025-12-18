##
## patient 3348 unknown
## patient 3354 nog aanvullen
pneumonia <- filter(EB_1.0, EB_1.0$group == "Ward")
pneumonia <- cohort_eb[cohort_eb$EB_id %in% pneumonia$ID, ]
table(pneumonia$Unknown)

##
one <- pneumonia$outcome_pat_cause1
two <- pneumonia$outcome_pat_cause2
three <- pneumonia$outcome_pat_cause3
all <- c(one, two, three)
print(table(all), quote = T)


##
table(pneumonia$outcome_pat_cause1)
pneumonia$strep <- ifelse(pneumonia$outcome_pat_cause1 == "Streptococcus pneumoniae", 1, 
                          ifelse(pneumonia$outcome_pat_cause2 == "Streptococcus pneumoniae", 1, 
                                 ifelse(pneumonia$outcome_pat_cause3 == "Streptococcus pneumoniae", 1, 0)))
pneumonia$strep <- ifelse(is.na(pneumonia$strep), 0, pneumonia$strep)

pneumonia$covid <- ifelse(pneumonia$outcome_pat_cause1 == "COVID-19", 1, 
                          ifelse(pneumonia$outcome_pat_cause2 == "COVID-19", 1, 
                                 ifelse(pneumonia$outcome_pat_cause3 == "COVID-19", 1, 0)))
pneumonia$covid <- ifelse(is.na(pneumonia$covid), 0, pneumonia$covid)

pneumonia$aureus <- ifelse(pneumonia$outcome_pat_cause1 == "Staphylococcus aureus", 1, 
                          ifelse(pneumonia$outcome_pat_cause2 == "Staphylococcus aureus", 1, 
                                 ifelse(pneumonia$outcome_pat_cause3 == "Staphylococcus aureus", 1, 0)))
pneumonia$aureus <- ifelse(is.na(pneumonia$aureus), 0, pneumonia$aureus)

pneumonia$influenza <- ifelse(pneumonia$outcome_pat_cause1 == "Influenza A virus", 1, 
                           ifelse(pneumonia$outcome_pat_cause2 == "Influenza A virus", 1, 
                                  ifelse(pneumonia$outcome_pat_cause3 == "Influenza A virus", 1, 0)))
pneumonia$influenza <- ifelse(is.na(pneumonia$influenza), 0, pneumonia$influenza)

pneumonia$h_influenz <- ifelse(pneumonia$outcome_pat_cause1 == "Haemophilus influenzae", 1, 
                              ifelse(pneumonia$outcome_pat_cause2 == "Haemophilus influenzae", 1, 
                                     ifelse(pneumonia$outcome_pat_cause3 == "Haemophilus influenzae", 1, 0)))
pneumonia$h_influenz <- ifelse(is.na(pneumonia$h_influenz), 0, pneumonia$h_influenz)

pneumonia$no_pathogen <- ifelse(pneumonia$outcome_pat_cause1 == "No causative agent found", 1, 0)

##
pneumonia$mixed_infections <- ifelse(pneumonia$outcome_pat_cause2 == 0, 0, 1)

## 
pneumonia$other_pathogen <- ifelse(pneumonia$no_pathogen == 1, 0,
                                   ifelse(pneumonia$strep == 1, 0, 
                                          ifelse(pneumonia$covid == 1, 0,
                                                 ifelse(pneumonia$influenza == 1, 0,
                                                        ifelse(pneumonia$h_influenz == 1, 0,
                                                               ifelse(pneumonia$aureus == 1, 0, 1))))))

##
pneumonia$strep <- ifelse(pneumonia$mixed_infections == 1, 0, pneumonia$strep)
pneumonia$influenza <- ifelse(pneumonia$mixed_infections == 1, 0, pneumonia$influenza)
pneumonia$covid <- ifelse(pneumonia$mixed_infections == 1, 0, pneumonia$covid)
pneumonia$aureus <- ifelse(pneumonia$mixed_infections == 1, 0, pneumonia$aureus)
pneumonia$h_influenz <- ifelse(pneumonia$mixed_infections == 1, 0, pneumonia$h_influenz)
pneumonia$other_pathogen <- ifelse(pneumonia$mixed_infections == 1, 0, pneumonia$other_pathogen)

##
setwd("//amc.intra/users/E/ehmichels/home/Desktop/CAP ICU no ICU/CAP ICU no ICU")
other_pathogens <- filter(pneumonia, pneumonia$other_pathogen == 1)
other_pathogens <- select(other_pathogens, "EB_id",  "outcome_pat_cause1", "outcome_def_pathogen1")
write.csv(other_pathogens, "other_pathogens.csv")

##
mixed_infection <- filter(pneumonia, pneumonia$mixed_infections == 1)
mixed_infection_short <- select(mixed_infection, "EB_id", "outcome_pat_cause1", "outcome_pat_cause2", "outcome_def_pathogen2", "outcome_pat_cause3", 
                                "outcome_def_pathogen3")
write.csv(mixed_infection_short, "mixed_infections.csv")

##
pneumonia$pathogen <- ifelse(pneumonia$mixed_infections == 1, "mixed_infection", 
                             ifelse(pneumonia$no_pathogen == 1, "No_pathogen", 
                             ifelse(pneumonia$strep == 1, "Strep", 
                                    ifelse(pneumonia$aureus == 1, "aureus", 
                                           ifelse(pneumonia$influenza == 1, "influenza", 
                                                  ifelse(pneumonia$h_influenz == 1, "haemo", 
                                                         ifelse(pneumonia$covid == 1, "covid", "other")))))))
table(pneumonia$pathogen)

## 65 have a pathogen
pneumonia_pathogen <- filter(pneumonia, !pneumonia$outcome_pat_cause1 == "No causative agent found")
pneumonia_pathogen

## binnen 48h (check)
## missing admission dates

## change to proper times
## missing dates checked => now okay
pneumonia_pathogen$base_date_admission <- pneumonia_pathogen$admission_dt1
pneumonia_pathogen$`outcome_pat_microbiology_results_1_Microbiological_test_1_Sample_date` <- as.Date(pneumonia_pathogen$`outcome_pat_microbiology_results_1_Microbiological_test_1_Sample_date`, format = "%d-%m-%Y")
pneumonia_pathogen$`outcome_pat_microbiology_results_1_Microbiological_test_2_Sample_date` <- as.Date(pneumonia_pathogen$`outcome_pat_microbiology_results_1_Microbiological_test_2_Sample_date`, format = "%d-%m-%Y")
pneumonia_pathogen$`outcome_pat_microbiology_results_1_Microbiological_test_3_Sample_date` <- as.Date(pneumonia_pathogen$`outcome_pat_microbiology_results_1_Microbiological_test_3_Sample_date`, format = "%d-%m-%Y")
pneumonia_pathogen$`outcome_pat_microbiology_results_1_Microbiological_test_4_Sample_date` <- as.Date(pneumonia_pathogen$`outcome_pat_microbiology_results_1_Microbiological_test_4_Sample_date`, format = "%d-%m-%Y")
pneumonia_pathogen$`outcome_pat_microbiology_results_1_Microbiological_test_5_Sample_date` <- as.Date(pneumonia_pathogen$`outcome_pat_microbiology_results_1_Microbiological_test_5_Sample_date`, format = "%d-%m-%Y")
pneumonia_pathogen$`outcome_pat_microbiology_results_1_Microbiological_test_6_Sample_date` <- as.Date(pneumonia_pathogen$`outcome_pat_microbiology_results_1_Microbiological_test_6_Sample_date`, format = "%d-%m-%Y")
pneumonia_pathogen$`outcome_pat_microbiology_results_1_Microbiological_test_7_Sample_date` <- as.Date(pneumonia_pathogen$`outcome_pat_microbiology_results_1_Microbiological_test_7_Sample_date`, format = "%d-%m-%Y")

## calculate time between admission and sampling
pneumonia_pathogen$admission_to_sample1 <- difftime(pneumonia_pathogen$`outcome_pat_microbiology_results_1_Microbiological_test_1_Sample_date`, pneumonia_pathogen$base_date_admission, units = "days")
pneumonia_pathogen$admission_to_sample2 <- difftime(pneumonia_pathogen$`outcome_pat_microbiology_results_1_Microbiological_test_2_Sample_date`, pneumonia_pathogen$base_date_admission, units = "days")
pneumonia_pathogen$admission_to_sample3 <- difftime(pneumonia_pathogen$`outcome_pat_microbiology_results_1_Microbiological_test_3_Sample_date`, pneumonia_pathogen$base_date_admission, units = "days")
pneumonia_pathogen$admission_to_sample4 <- difftime(pneumonia_pathogen$`outcome_pat_microbiology_results_1_Microbiological_test_4_Sample_date`, pneumonia_pathogen$base_date_admission, units = "days")
pneumonia_pathogen$admission_to_sample5 <- difftime(pneumonia_pathogen$`outcome_pat_microbiology_results_1_Microbiological_test_5_Sample_date`, pneumonia_pathogen$base_date_admission, units = "days")
pneumonia_pathogen$admission_to_sample6 <- difftime(pneumonia_pathogen$`outcome_pat_microbiology_results_1_Microbiological_test_6_Sample_date`, pneumonia_pathogen$base_date_admission, units = "days")
pneumonia_pathogen$admission_to_sample7 <- difftime(pneumonia_pathogen$`outcome_pat_microbiology_results_1_Microbiological_test_7_Sample_date`, pneumonia_pathogen$base_date_admission, units = "days")


pneumonia_pathogen_1 <- select(pneumonia_pathogen, "EB_id", "outcome_pat_cause1", "outcome_pat_cause2", "outcome_pat_cause3", "admission_to_sample1",
                    "base_date_admission", "outcome_pat_microbiology_results_1_Microbiological_test_1_Sample_origin",
                    "outcome_pat_microbiology_results_1_Microbiological_test_1_Result", "outcome_pat_microbiology_results_1_Microbiological_test_1_BACTERIA_found",
                    "outcome_pat_microbiology_results_1_Microbiological_test_1_VIRUS_found", "outcome_pat_microbiology_results_1_Microbiological_test_1_FUNGI_found", 
                    "outcome_pat_microbiology_results_1_Microbiological_test_1_Specify_other_pathogen")

pneumonia_pathogen_2 <- select(pneumonia_pathogen, "EB_id", "outcome_pat_cause1", "outcome_pat_cause2", "outcome_pat_cause3", "admission_to_sample2",
                      "base_date_admission","outcome_pat_microbiology_results_1_Microbiological_test_2_Sample_origin",
                      "outcome_pat_microbiology_results_1_Microbiological_test_2_Result", "outcome_pat_microbiology_results_1_Microbiological_test_2_BACTERIA_found",
                      "outcome_pat_microbiology_results_1_Microbiological_test_2_VIRUS_found", "outcome_pat_microbiology_results_1_Microbiological_test_2_FUNGI_found", 
                      "outcome_pat_microbiology_results_1_Microbiological_test_2_Specify_other_pathogen")

pneumonia_pathogen_3 <- select(pneumonia_pathogen, "EB_id", "outcome_pat_cause1", "outcome_pat_cause2", "outcome_pat_cause3", "admission_to_sample3",
                      "base_date_admission","outcome_pat_microbiology_results_1_Microbiological_test_3_Sample_origin",
                      "outcome_pat_microbiology_results_1_Microbiological_test_3_Result", "outcome_pat_microbiology_results_1_Microbiological_test_3_BACTERIA_found",
                      "outcome_pat_microbiology_results_1_Microbiological_test_3_VIRUS_found", "outcome_pat_microbiology_results_1_Microbiological_test_3_FUNGI_found", 
                      "outcome_pat_microbiology_results_1_Microbiological_test_3_Specify_other_pathogen")

pneumonia_pathogen_4 <- select(pneumonia_pathogen, "EB_id", "outcome_pat_cause1", "outcome_pat_cause2", "outcome_pat_cause3", "admission_to_sample4",
                      "base_date_admission","outcome_pat_microbiology_results_1_Microbiological_test_4_Sample_origin",
                      "outcome_pat_microbiology_results_1_Microbiological_test_4_Result", "outcome_pat_microbiology_results_1_Microbiological_test_4_BACTERIA_found",
                      "outcome_pat_microbiology_results_1_Microbiological_test_4_VIRUS_found", "outcome_pat_microbiology_results_1_Microbiological_test_4_FUNGI_found", 
                      "outcome_pat_microbiology_results_1_Microbiological_test_4_Specify_other_pathogen")

pneumonia_pathogen_5 <- select(pneumonia_pathogen, "EB_id", "outcome_pat_cause1", "outcome_pat_cause2", "outcome_pat_cause3", "admission_to_sample5",
                      "base_date_admission","outcome_pat_microbiology_results_1_Microbiological_test_5_Sample_origin",
                      "outcome_pat_microbiology_results_1_Microbiological_test_5_Result", "outcome_pat_microbiology_results_1_Microbiological_test_5_BACTERIA_found",
                      "outcome_pat_microbiology_results_1_Microbiological_test_5_VIRUS_found", "outcome_pat_microbiology_results_1_Microbiological_test_5_FUNGI_found", 
                      "outcome_pat_microbiology_results_1_Microbiological_test_5_Specify_other_pathogen")

pneumonia_pathogen_6 <- select(pneumonia_pathogen, "EB_id", "outcome_pat_cause1", "outcome_pat_cause2", "outcome_pat_cause3", "admission_to_sample6",
                      "base_date_admission","outcome_pat_microbiology_results_1_Microbiological_test_6_Sample_origin",
                      "outcome_pat_microbiology_results_1_Microbiological_test_6_Result", "outcome_pat_microbiology_results_1_Microbiological_test_6_BACTERIA_found",
                      "outcome_pat_microbiology_results_1_Microbiological_test_6_VIRUS_found", "outcome_pat_microbiology_results_1_Microbiological_test_6_FUNGI_found", 
                      "outcome_pat_microbiology_results_1_Microbiological_test_6_Specify_other_pathogen")

pneumonia_pathogen_7 <- select(pneumonia_pathogen, "EB_id", "outcome_pat_cause1", "outcome_pat_cause2", "outcome_pat_cause3", "admission_to_sample7",
                      "base_date_admission", "outcome_pat_microbiology_results_1_Microbiological_test_7_Sample_origin",
                      "outcome_pat_microbiology_results_1_Microbiological_test_7_Result", "outcome_pat_microbiology_results_1_Microbiological_test_7_BACTERIA_found",
                      "outcome_pat_microbiology_results_1_Microbiological_test_7_VIRUS_found", "outcome_pat_microbiology_results_1_Microbiological_test_7_FUNGI_found", 
                      "outcome_pat_microbiology_results_1_Microbiological_test_7_Specify_other_pathogen")

##
x1 <- pneumonia_pathogen_1[is.na(pneumonia_pathogen_1$admission_to_sample1), ]
x2 <- pneumonia_pathogen_2[is.na(pneumonia_pathogen_2$admission_to_sample2), ]
x3 <- pneumonia_pathogen_3[is.na(pneumonia_pathogen_3$admission_to_sample3), ]
x4 <- pneumonia_pathogen_4[is.na(pneumonia_pathogen_4$admission_to_sample4), ]
x5 <- pneumonia_pathogen_5[is.na(pneumonia_pathogen_5$admission_to_sample5), ]
x6 <- pneumonia_pathogen_6[is.na(pneumonia_pathogen_6$admission_to_sample6), ]
x7 <- pneumonia_pathogen_7[is.na(pneumonia_pathogen_7$admission_to_sample7), ]


## 8 missing result + more missing Sample_date
## check al to high Sample_date
check<- filter(pneumonia_pathogen_1, is.na(pneumonia_pathogen_1$`outcome_pat_microbiology_results_1_Microbiological_test_1_Result`))
pneumonia_pathogen_1 <- filter(pneumonia_pathogen_1, 
                                !pneumonia_pathogen_1$`outcome_pat_microbiology_results_1_Microbiological_test_1_Result` == "No growth/pathogen found")
pneumonia_pathogen_admission_1 <- filter(pneumonia_pathogen_1, is.na(pneumonia_pathogen_1$admission_to_sample1) |  pneumonia_pathogen_1$admission_to_sample1 >= 2)

colnames(pneumonia_pathogen_admission_1)[colnames(pneumonia_pathogen_admission_1) == "admission_to_sample1"] <- "admission_to_sample"
colnames(pneumonia_pathogen_admission_1)[colnames(pneumonia_pathogen_admission_1) == "outcome_pat_microbiology_results_1_Microbiological_test_1_Result"] <- "Result"
colnames(pneumonia_pathogen_admission_1)[colnames(pneumonia_pathogen_admission_1) == "outcome_pat_microbiology_results_1_Microbiological_test_1_BACTERIA_found"] <- "Bacteria"
colnames(pneumonia_pathogen_admission_1)[colnames(pneumonia_pathogen_admission_1) == "outcome_pat_microbiology_results_1_Microbiological_test_1_VIRUS_found"] <- "Virus"
colnames(pneumonia_pathogen_admission_1)[colnames(pneumonia_pathogen_admission_1) == "outcome_pat_microbiology_results_1_Microbiological_test_1_FUNGI_found"] <- "Fungi"
colnames(pneumonia_pathogen_admission_1)[colnames(pneumonia_pathogen_admission_1) == "outcome_pat_microbiology_results_1_Microbiological_test_1_Specify_other_pathogen"] <- "Other"
colnames(pneumonia_pathogen_admission_1)[colnames(pneumonia_pathogen_admission_1) == "outcome_pat_microbiology_results_1_Microbiological_test_1_Sample_origin"] <- "Origin"


##
pneumonia_pathogen_2 <- filter(pneumonia_pathogen_2, 
                                 !pneumonia_pathogen_2$`outcome_pat_microbiology_results_1_Microbiological_test_2_Result` == "No growth/pathogen found")
pneumonia_pathogen_admission_2 <- filter(pneumonia_pathogen_2, is.na(pneumonia_pathogen_2$admission_to_sample2) |  pneumonia_pathogen_2$admission_to_sample2 >= 2)

colnames(pneumonia_pathogen_admission_2)[colnames(pneumonia_pathogen_admission_2) == "admission_to_sample2"] <- "admission_to_sample"
colnames(pneumonia_pathogen_admission_2)[colnames(pneumonia_pathogen_admission_2) == "outcome_pat_microbiology_results_1_Microbiological_test_2_Result"] <- "Result"
colnames(pneumonia_pathogen_admission_2)[colnames(pneumonia_pathogen_admission_2) == "outcome_pat_microbiology_results_1_Microbiological_test_2_BACTERIA_found"] <- "Bacteria"
colnames(pneumonia_pathogen_admission_2)[colnames(pneumonia_pathogen_admission_2) == "outcome_pat_microbiology_results_1_Microbiological_test_2_VIRUS_found"] <- "Virus"
colnames(pneumonia_pathogen_admission_2)[colnames(pneumonia_pathogen_admission_2) == "outcome_pat_microbiology_results_1_Microbiological_test_2_FUNGI_found"] <- "Fungi"
colnames(pneumonia_pathogen_admission_2)[colnames(pneumonia_pathogen_admission_2) == "outcome_pat_microbiology_results_1_Microbiological_test_2_Specify_other_pathogen"] <- "Other"
colnames(pneumonia_pathogen_admission_2)[colnames(pneumonia_pathogen_admission_2) == "outcome_pat_microbiology_results_1_Microbiological_test_2_Sample_origin"] <- "Origin"

##
pneumonia_pathogen_3 <- filter(pneumonia_pathogen_3, 
                                 !pneumonia_pathogen_3$`outcome_pat_microbiology_results_1_Microbiological_test_3_Result` == "No growth/pathogen found")
pneumonia_pathogen_admission_3 <- filter(pneumonia_pathogen_3, is.na(pneumonia_pathogen_3$admission_to_sample3) |  pneumonia_pathogen_3$admission_to_sample3 >= 2)

colnames(pneumonia_pathogen_admission_3)[colnames(pneumonia_pathogen_admission_3) == "admission_to_sample3"] <- "admission_to_sample"
colnames(pneumonia_pathogen_admission_3)[colnames(pneumonia_pathogen_admission_3) == "outcome_pat_microbiology_results_1_Microbiological_test_3_Result"] <- "Result"
colnames(pneumonia_pathogen_admission_3)[colnames(pneumonia_pathogen_admission_3) == "outcome_pat_microbiology_results_1_Microbiological_test_3_BACTERIA_found"] <- "Bacteria"
colnames(pneumonia_pathogen_admission_3)[colnames(pneumonia_pathogen_admission_3) == "outcome_pat_microbiology_results_1_Microbiological_test_3_VIRUS_found"] <- "Virus"
colnames(pneumonia_pathogen_admission_3)[colnames(pneumonia_pathogen_admission_3) == "outcome_pat_microbiology_results_1_Microbiological_test_3_FUNGI_found"] <- "Fungi"
colnames(pneumonia_pathogen_admission_3)[colnames(pneumonia_pathogen_admission_3) == "outcome_pat_microbiology_results_1_Microbiological_test_3_Specify_other_pathogen"] <- "Other"
colnames(pneumonia_pathogen_admission_3)[colnames(pneumonia_pathogen_admission_3) == "outcome_pat_microbiology_results_1_Microbiological_test_3_Sample_origin"] <- "Origin"


##
pneumonia_pathogen_4 <- filter(pneumonia_pathogen_4, 
                                 !pneumonia_pathogen_4$`outcome_pat_microbiology_results_1_Microbiological_test_4_Result` == "No growth/pathogen found")
pneumonia_pathogen_admission_4 <- filter(pneumonia_pathogen_4, is.na(pneumonia_pathogen_4$admission_to_sample4) |  pneumonia_pathogen_4$admission_to_sample4 >= 2)

colnames(pneumonia_pathogen_admission_4)[colnames(pneumonia_pathogen_admission_4) == "admission_to_sample4"] <- "admission_to_sample"
colnames(pneumonia_pathogen_admission_4)[colnames(pneumonia_pathogen_admission_4) == "outcome_pat_microbiology_results_1_Microbiological_test_4_Result"] <- "Result"
colnames(pneumonia_pathogen_admission_4)[colnames(pneumonia_pathogen_admission_4) == "outcome_pat_microbiology_results_1_Microbiological_test_4_BACTERIA_found"] <- "Bacteria"
colnames(pneumonia_pathogen_admission_4)[colnames(pneumonia_pathogen_admission_4) == "outcome_pat_microbiology_results_1_Microbiological_test_4_VIRUS_found"] <- "Virus"
colnames(pneumonia_pathogen_admission_4)[colnames(pneumonia_pathogen_admission_4) == "outcome_pat_microbiology_results_1_Microbiological_test_4_FUNGI_found"] <- "Fungi"
colnames(pneumonia_pathogen_admission_4)[colnames(pneumonia_pathogen_admission_4) == "outcome_pat_microbiology_results_1_Microbiological_test_4_Specify_other_pathogen"] <- "Other"
colnames(pneumonia_pathogen_admission_4)[colnames(pneumonia_pathogen_admission_4) == "outcome_pat_microbiology_results_1_Microbiological_test_4_Sample_origin"] <- "Origin"


##
pneumonia_pathogen_5 <- filter(pneumonia_pathogen_5, 
                                 !pneumonia_pathogen_5$`outcome_pat_microbiology_results_1_Microbiological_test_5_Result` == "No growth/pathogen found")
pneumonia_pathogen_admission_5 <- filter(pneumonia_pathogen_5, is.na(pneumonia_pathogen_5$admission_to_sample5) |  pneumonia_pathogen_5$admission_to_sample5 >= 2)

colnames(pneumonia_pathogen_admission_5)[colnames(pneumonia_pathogen_admission_5) == "admission_to_sample5"] <- "admission_to_sample"
colnames(pneumonia_pathogen_admission_5)[colnames(pneumonia_pathogen_admission_5) == "outcome_pat_microbiology_results_1_Microbiological_test_5_Result"] <- "Result"
colnames(pneumonia_pathogen_admission_5)[colnames(pneumonia_pathogen_admission_5) == "outcome_pat_microbiology_results_1_Microbiological_test_5_BACTERIA_found"] <- "Bacteria"
colnames(pneumonia_pathogen_admission_5)[colnames(pneumonia_pathogen_admission_5) == "outcome_pat_microbiology_results_1_Microbiological_test_5_VIRUS_found"] <- "Virus"
colnames(pneumonia_pathogen_admission_5)[colnames(pneumonia_pathogen_admission_5) == "outcome_pat_microbiology_results_1_Microbiological_test_5_FUNGI_found"] <- "Fungi"
colnames(pneumonia_pathogen_admission_5)[colnames(pneumonia_pathogen_admission_5) == "outcome_pat_microbiology_results_1_Microbiological_test_5_Specify_other_pathogen"] <- "Other"
colnames(pneumonia_pathogen_admission_5)[colnames(pneumonia_pathogen_admission_5) == "outcome_pat_microbiology_results_1_Microbiological_test_5_Sample_origin"] <- "Origin"


##
pneumonia_pathogen_6 <- filter(pneumonia_pathogen_6, 
                                 !pneumonia_pathogen_6$`outcome_pat_microbiology_results_1_Microbiological_test_6_Result` == "No growth/pathogen found")
pneumonia_pathogen_admission_6 <- filter(pneumonia_pathogen_6, is.na(pneumonia_pathogen_6$admission_to_sample6) |  pneumonia_pathogen_6$admission_to_sample6 >= 2)

colnames(pneumonia_pathogen_admission_6)[colnames(pneumonia_pathogen_admission_6) == "admission_to_sample6"] <- "admission_to_sample"
colnames(pneumonia_pathogen_admission_6)[colnames(pneumonia_pathogen_admission_6) == "outcome_pat_microbiology_results_1_Microbiological_test_6_Result"] <- "Result"
colnames(pneumonia_pathogen_admission_6)[colnames(pneumonia_pathogen_admission_6) == "outcome_pat_microbiology_results_1_Microbiological_test_6_BACTERIA_found"] <- "Bacteria"
colnames(pneumonia_pathogen_admission_6)[colnames(pneumonia_pathogen_admission_6) == "outcome_pat_microbiology_results_1_Microbiological_test_6_VIRUS_found"] <- "Virus"
colnames(pneumonia_pathogen_admission_6)[colnames(pneumonia_pathogen_admission_6) == "outcome_pat_microbiology_results_1_Microbiological_test_6_FUNGI_found"] <- "Fungi"
colnames(pneumonia_pathogen_admission_6)[colnames(pneumonia_pathogen_admission_6) == "outcome_pat_microbiology_results_1_Microbiological_test_6_Specify_other_pathogen"] <- "Other"
colnames(pneumonia_pathogen_admission_6)[colnames(pneumonia_pathogen_admission_6) == "outcome_pat_microbiology_results_1_Microbiological_test_6_Sample_origin"] <- "Origin"


##
pneumonia_pathogen_7 <- filter(pneumonia_pathogen_7,
                                 !pneumonia_pathogen_7$`outcome_pat_microbiology_results_1_Microbiological_test_7_Result` == "No growth/pathogen found")
pneumonia_pathogen_admission_7 <- filter(pneumonia_pathogen_7, is.na(pneumonia_pathogen_7$admission_to_sample7) |  pneumonia_pathogen_7$admission_to_sample7 >= 2)

colnames(pneumonia_pathogen_admission_7)[colnames(pneumonia_pathogen_admission_7) == "admission_to_sample7"] <- "admission_to_sample"
colnames(pneumonia_pathogen_admission_7)[colnames(pneumonia_pathogen_admission_7) == "outcome_pat_microbiology_results_1_Microbiological_test_7_Result"] <- "Result"
colnames(pneumonia_pathogen_admission_7)[colnames(pneumonia_pathogen_admission_7) == "outcome_pat_microbiology_results_1_Microbiological_test_7_BACTERIA_found"] <- "Bacteria"
colnames(pneumonia_pathogen_admission_7)[colnames(pneumonia_pathogen_admission_7) == "outcome_pat_microbiology_results_1_Microbiological_test_7_VIRUS_found"] <- "Virus"
colnames(pneumonia_pathogen_admission_7)[colnames(pneumonia_pathogen_admission_7) == "outcome_pat_microbiology_results_1_Microbiological_test_7_FUNGI_found"] <- "Fungi"
colnames(pneumonia_pathogen_admission_7)[colnames(pneumonia_pathogen_admission_7) == "outcome_pat_microbiology_results_1_Microbiological_test_7_Specify_other_pathogen"] <- "Other"
colnames(pneumonia_pathogen_admission_7)[colnames(pneumonia_pathogen_admission_7) == "outcome_pat_microbiology_results_1_Microbiological_test_7_Sample_origin"] <- "Origin"

##

## 
all_pathogen <- rbind(pneumonia_pathogen_admission_1, pneumonia_pathogen_admission_2, pneumonia_pathogen_admission_3, 
                      pneumonia_pathogen_admission_4, pneumonia_pathogen_admission_5, pneumonia_pathogen_admission_6,pneumonia_pathogen_admission_7)


##
writexl::write_xlsx(all_pathogen, "all_acquired_pathogen.xlsx")
suspected_pathogen <- filter(all_pathogen, all_pathogen$Result == "Suspected pathogen found")

#
tiwce <- twice[order(twice$`EB_id`), ]
setwd("/Users/erikmichels/Desktop/PhD/Stimulaties")
writexl::write_xlsx(once, "one_pathogen.xlsx")
writexl::write_xlsx(twice, "twice_pathogen.xlsx")
writexl::write_xlsx(all_pathogen, "all_pathogen.xlsx")

## only suspected 
## manually inspected "Atypical" have no pathgoen
suspected_pathogen <- filter(all_pathogen, all_pathogen$Result == "Suspected pathogen found")
suspected_pathogen <- suspected_pathogen[order(suspected_pathogen$`EB_id`), ]
writexl::write_xlsx(suspected_pathogen, "suspected_pathogen.xlsx")
 
#
## only 3 patients