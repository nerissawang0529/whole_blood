rm(list = ls())

studydata_patients <- read.csv("Original_data/studydata_patients.csv")
#studydata_HV <- read.csv("Original_data/studydata_HV.csv")
#studydata_combined <- rbind(studydata_patients, studydata_HV)

Luminex_CAP_COVID_HV <- read.csv("Original_data/Luminex_CAP_COVID_HV.csv")
merged_data <- merge(Luminex_CAP_COVID_HV, studydata_patients, 
                     by.x = "ID", by.y = "EB_id")



mdata <- merged_data[,1:12]
mdata$Day <- NULL

library(reshape2)
ldata <- reshape2::melt(mdata, id=c("ID","stimulation") )
ldata$level <- log10(ldata$value)

hist(ldata$level)
library(ggplot2)
ggplot(ldata, aes(x=level, fill=stimulation))+
  geom_histogram(alpha=0.8, binwidth = 0.2)+
  facet_wrap(~variable)+
  theme_classic()

library(lme4)
library(lmerTest)

ldata$stimulation <- as.factor(ldata$stimulation)
ldata$stimulation <- factor(ldata$stimulation, levels = c("M","KP","LPS","PP"))

ldata_keep <- ldata

ldata <- ldata[!is.na(ldata$level),]


mod0 <- lmerTest::lmer(level ~ variable + stimulation + (1|ID) , data=ldata)

summary(mod0)


mod <- lmerTest::lmer(level ~ variable * stimulation + (1|ID) , data=ldata)

summary(mod)


anova(mod0, mod)


residualdf <- residuals(mod)

predictdf <- predict(mod)


predictdf_fixedonly <- predict(mod, re.form=NA, newdata=ldata_keep )
residualdf_fixedonly <- data.frame(ldata_keep, predicted=predictdf_fixedonly)

residualdf_fixedonly[residualdf_fixedonly$ID==3005 & 
                       residualdf_fixedonly$variable=="IL_1beta" &
                       residualdf_fixedonly$stimulation=="PP",]

residualdf_fixedonly[residualdf_fixedonly$ID==3005 & 
                       residualdf_fixedonly$variable=="IL_1beta" &
                       residualdf_fixedonly$stimulation=="PP",]

residualdf_fixedonly[is.na(residualdf_fixedonly$value),]




residualdf <- data.frame(ldata, predicted=predict(mod),  resid= residuals(mod))

ggplot(residualdf, aes(x=resid, fill=stimulation))+
  geom_histogram(alpha=0.8)+
  facet_wrap(~variable)+
  theme_classic()


#### make wide residuals 

residualdf0 <- residualdf
residualdf0$level <- NULL
residualdf0$predicted <- NULL


resid_w <- reshape2::dcast(residualdf0, ID ~ stimulation + variable, value.var="resid" )

row.names(resid_w) <- resid_w$ID
resid_w$ID <- NULL


resid_w["3005", "PP_IL_1beta"] <- 0
resid_w["3350", "PP_CCL3"] <- 0


pca_res <- prcomp(resid_w)

library(ggfortify)

autoplot(pca_res)


autoplot(pca_res, 
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3)


kmeansres <- kmeans(pca_res$x[,1:2], centers = 2) 

pat_clust <- data.frame(cluster=as.factor(kmeansres$cluster))

autoplot(pca_res, data=pat_clust, colour="cluster",
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3)

autoplot(pca_res, data=pat_clust, colour="cluster")
table(pat_clust$cluster)
pat_clust$ID <- rownames(pat_clust)

##table1
studydata_patients <- read.csv("Original_data/studydata_patients.csv")
studydata_patients_0 <- merge(studydata_patients,pat_clust,by.x = "EB_id", by.y = "ID")

## table
allvars <- c("inclusion_hospital", "sampling_time",		"age_yrs",	"gender",	"ethnic_group#White",	
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
             "hospdeath",	"mortality_d30", "mortality_d90","pathogen","Gram_negative","Gram_positive","Virus","Fungi","cluster")

catvars <- c("inclusion_hospital", "gender",	"ethnic_group#White",	
             "hypertension",	"cpd",	"COPD",	"diabetes",	
             "ccd","ckd",	"mneoplasm",	"immune_sup",	"cnd",	
             "Oxygen_therapy_1", 
             "antibiotic_seven_days", 
             "ICU_Medium_Care_admission_1",	 
             "Non_invasive_ventilation_1",	"Invasive_ventilation_1",	
             "bacteremia",
             "hospdeath",	"mortality_d30", "mortality_d90", "pathogen_cultured","pathogen","Gram_negative","Gram_positive","Virus","Fungi")

nonnormal <- c("sampling_time",	"age_yrs", "BMI", "symptom_days",  "CCI", "qSOFA_score",	"MEWS_score",	"CURB_score",	"PSI_new",
               "crp_1_1",	 "WBC_2_1",	"Neutrophil_unit_1",	"Lymphocyte_1_1",	
               "Creatinine_value_1",
               "Platelets_value_1",	
               "Blood_Urea_Nitrogen_value_1",	"length_of_oxygen","length_of_stay" ,"icu_stay", "lenght_of_intubation")

tab2     <- CreateTableOne(
  vars        = allvars, 
  data        = studydata_patients_0, 
  strata = "cluster",
  factorVars  = catvars,
  test        = TRUE)

print(tab2, nonnormal = nonnormal, quote = TRUE, noSpaces = TRUE, smd = T, missing = T)



> dim(TNF_M)
[1] 221   5
> TNF_L <- ldata[ldata$variable=="TNF" & ldata$stimulation=="LPS",]
> dim(TNF_L)
[1] 221   5
> hist(TNF_L)
> hist(TNF_L$level)
> hist(TNF_M$level)
> IL6_M <- ldata[ldata$variable=="IL_6" & ldata$stimulation=="M",]
> hist(IL6_M$level)
> diptest::dip(IL6_M$level)
[1] 0.03088243
> diptest::dip(IL8_M$level)
[1] 0.03528679
> diptest::dip(TNF_M$level)
[1] 0.01334001
> diptest::dip(IL_1RA_M$level)
[1] 0.01304976
> diptest::dip(CCL2_M$level)
[1] 0.02199482
> diptest::dip(CCL3_M$level)
[1] 0.01900182
> diptest::dip(IL_1beta_M$level)
[1] 0.03010187
> diptest::dip(IL_10_M$level)
[1] 0.02415774