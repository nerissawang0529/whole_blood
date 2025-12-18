rm(list = ls())
## install.packages("pacman")
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
               tableone, kmed, survival, survminer, skimr, lattice, reshape2, grDevices, ggrepel)
library("variancePartition")

studydata_patients <- read.csv("Original_data/studydata_patients.csv")
studydata_patients <- studydata_patients[!studydata_patients$EB_id %in% c(3348, 1167, 1237, 1238, 3273, 3292, 3321), ]
#the MEWS is NA

Luminex_CAP_COVID_HV <- read.csv("Original_data/Luminex_CAP_COVID_HV_log.csv")
Luminex_CAP_COVID_HV <- Luminex_CAP_COVID_HV[!Luminex_CAP_COVID_HV$ID %in% c(3348, 1167, 1237, 1238, 3273, 3292, 3321), ]
merged_data <- merge(Luminex_CAP_COVID_HV, studydata_patients, 
                     by.x = "ID", by.y = "EB_id")

#ready for exp2####
mdata <- merged_data[,c(1,12,14:22)]
mdata <- mdata %>% filter(!stimulation == "M")
mdata <- na.omit(mdata)

library(reshape2)
ldata <- reshape2::melt(mdata, id=c("ID","stimulation") )

# Convert ldata to a data.table
setDT(ldata)
library(reshape2)

exp2 <- dcast(ldata, variable ~ ID + stimulation, value.var = "value")
exp2 <- as.data.frame(exp2)
exp2$variable <- as.character(exp2$variable)
rownames(exp2) <- exp2$variable
exp2 <- exp2[, -1]  #####

# bsinfo####
library(dplyr)
bsinfo <- merged_data[,c("ID", "stimulation", "age_yrs", "gender", "MEWS_score")]
colnames(bsinfo) <- c("ID", "Stimulation", "Age", "Sex", "MEWS")
bsinfo <- bsinfo %>%
  distinct() %>%
  mutate(Stimulation = factor(Stimulation),
         Sex = factor(Sex),
         ID = factor(ID))  

# Ensure rownames match exp2 column names: "ID_Stimulation"
rownames(bsinfo) <- paste(bsinfo$ID, bsinfo$Stimulation, sep = "_")

# Convert MEWS to numeric if needed
bsinfo$MEWS <- as.numeric(bsinfo$MEWS)
colnames(bsinfo)[colnames(bsinfo) == "ID"] <- "Individual"

# Check structure
str(bsinfo)

#form####
form <- ~ (1 | Stimulation) + Age + (1 | Sex) + MEWS + (1 | Individual)
# Reorder bsinfo to match the column order in exp2
bsinfo <- bsinfo[match(colnames(exp2), rownames(bsinfo)), ]


varPart <- fitExtractVarPartModel(exp2, form, bsinfo)

plotVarPart(varPart)

vp <- sortCols(varPart)

plotPercentBars(vp)

vp
