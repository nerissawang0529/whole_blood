rm(list = ls())
## install.packages("pacman")
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
               tableone, kmed, survival, survminer, skimr, lattice, reshape2, grDevices, ggrepel)
library("variancePartition")

studydata_patients <- read.csv("Original_data/studydata_patients.csv")
studydata_patients <- studydata_patients[!studydata_patients$EB_id %in% c(3348, 1167, 1237, 1238, 3273, 3292, 3321), ]
#the MEWS is NA

Luminex_CAP_COVID_HV <- read.csv("Original_data/Luminex_CAP_COVID_HV.csv")
Luminex_CAP_COVID_HV <- Luminex_CAP_COVID_HV[!Luminex_CAP_COVID_HV$ID %in% c(3348, 1167, 1237, 1238, 3273, 3292, 3321), ]
merged_data <- merge(Luminex_CAP_COVID_HV, studydata_patients, 
                     by.x = "ID", by.y = "EB_id")

#ready for exp2####
mdata <- merged_data[,1:12]
mdata <- na.omit(mdata)


mdata$Day <- NULL

library(reshape2)
ldata <- reshape2::melt(mdata, id=c("ID","stimulation") )
ldata$level <- log10(ldata$value)

ldata$value <- NULL
# Convert ldata to a data.table
setDT(ldata)
library(reshape2)

exp2 <- dcast(ldata, variable ~ ID + stimulation, value.var = "level")
exp2 <- as.data.frame(exp2)
exp2$variable <- as.character(exp2$variable)
rownames(exp2) <- exp2$variable
exp2 <- exp2[, -1]  #####

library(dplyr)


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


#### 
# Convert a categorical variable into indicator variables (make the stimulation into 4, the result is not reliable)


bsinfo$Stimulationf <- as.factor(bsinfo$Stimulation)
bsinfo$Stimulationf <- factor(bsinfo$Stimulationf, levels=c("M","LPS","PP","KP"))
dummy_matrix <- model.matrix(~ Stimulationf - 1, data = bsinfo)



# Convert matrix to dataframe for easier manipulation
dummy_df <- as.data.frame(dummy_matrix)
# Convert 0/1 to factor
dummy_df <- dummy_df %>%
  mutate(across(everything(), ~ as.factor(.)))
# View first few rows
head(dummy_df)

form2 <- ~ (1 | StimulationfLPS) + (1 | StimulationfPP) + (1 | StimulationfKP) + Age + (1 | Sex) + MEWS + (1 | Individual)

bsinfo2 <- cbind(bsinfo, dummy_df)

varPart2 <- fitExtractVarPartModel(exp2, form2, bsinfo2)

plotVarPart(varPart2)

vp2 <- sortCols(varPart2)

plotPercentBars(vp2)




##find the difference of didfferent groups
Luminex_CAP_COVID_HV <- read.csv("Original_data/Luminex_CAP_COVID_HV.csv")

studydata_patients <- read.csv("Original_data/studydata_patients.csv")
studydata_HV <- read.csv("Original_data/studydata_HV.csv")
studydata <- rbind(studydata_patients, studydata_HV)

merged_data <- merge(Luminex_CAP_COVID_HV, studydata, 
                     by.x = "ID", by.y = "EB_id")

table(merged_data$group)


#ready for exp2####
mdata <- merged_data[,1:12]
mdata <- na.omit(mdata)


mdata$Day <- NULL

library(reshape2)
ldata <- reshape2::melt(mdata, id=c("ID","stimulation") )
ldata$level <- log10(ldata$value)

table(merged_data$ID %in% ldata$ID)



ldata$value <- NULL
# Convert ldata to a data.table
setDT(ldata)
library(reshape2)

exp2 <- dcast(ldata, variable ~ ID + stimulation, value.var = "level")
exp2 <- as.data.frame(exp2)
exp2$variable <- as.character(exp2$variable)
rownames(exp2) <- exp2$variable
exp2 <- exp2[, -1]  #####


# bsinfo####
library(dplyr)
bsinfo <- merged_data[,c("ID", "stimulation", "age_yrs", "gender","group")]
colnames(bsinfo) <- c("ID", "Stimulation", "Age", "Sex","group")
bsinfo <- bsinfo %>%
  distinct() %>%
  mutate(Stimulation = factor(Stimulation),
         Sex = factor(Sex),
         ID = factor(ID),
         group = factor(group))  


rownames(bsinfo) <- paste(bsinfo$ID, bsinfo$Stimulation, sep = "_")

# Convert MEWS to numeric if needed
colnames(bsinfo)[colnames(bsinfo) == "ID"] <- "Individual"

colnames(bsinfo)[colnames(bsinfo) == "group"] <- "CAPvControl"
table(bsinfo$CAPvControl )


form <- ~ (1 | Stimulation) + Age + (1 | Sex) + (1|CAPvControl) + (1 | Individual)
# Reorder bsinfo to match the column order in exp2
bsinfo <- bsinfo[match(colnames(exp2), rownames(bsinfo)), ]


varPart <- fitExtractVarPartModel(exp2, form, bsinfo)

plotVarPart(varPart)

vp <- sortCols(varPart)

plotPercentBars(vp)

vp



### try for linear model 
Luminex_CAP_COVID_HV <- read.csv("Original_data/Luminex_CAP_COVID_HV.csv")

studydata_patients <- read.csv("Original_data/studydata_patients.csv")
studydata_HV <- read.csv("Original_data/studydata_HV.csv")
studydata <- rbind(studydata_patients, studydata_HV)

merged_data <- merge(Luminex_CAP_COVID_HV, studydata, 
                     by.x = "ID", by.y = "EB_id")

table(merged_data$group)


#ready for exp2####
mdata <- merged_data[,1:12]
mdata <- na.omit(mdata)


mdata <- mdata[mdata$stimulation !="M",]

mdata$Day <- NULL

library(reshape2)
ldata <- reshape2::melt(mdata, id=c("ID","stimulation") )
ldata$level <- log10(ldata$value)

table(merged_data$ID %in% ldata$ID)



ldata$value <- NULL
# Convert ldata to a data.table
setDT(ldata)
library(reshape2)

exp2 <- dcast(ldata, variable ~ ID + stimulation, value.var = "level")
exp2 <- as.data.frame(exp2)
exp2$variable <- as.character(exp2$variable)
rownames(exp2) <- exp2$variable
exp2 <- exp2[, -1]  #####


# bsinfo####
library(dplyr)
bsinfo <- merged_data[,c("ID", "stimulation", "age_yrs", "gender","group")]
colnames(bsinfo) <- c("ID", "Stimulation", "Age", "Sex","group")
bsinfo <- bsinfo %>%
  distinct() %>%
  mutate(Stimulation = factor(Stimulation),
         Sex = factor(Sex),
         ID = factor(ID),
         group = factor(group))  


rownames(bsinfo) <- paste(bsinfo$ID, bsinfo$Stimulation, sep = "_")

# Convert MEWS to numeric if needed
colnames(bsinfo)[colnames(bsinfo) == "ID"] <- "Individual"

colnames(bsinfo)[colnames(bsinfo) == "group"] <- "CAPvControl"
table(bsinfo$CAPvControl )


form <- ~ (1 | Stimulation) + Age + (1 | Sex) + (1|CAPvControl) 
# Reorder bsinfo to match the column order in exp2
bsinfo <- bsinfo[match(colnames(exp2), rownames(bsinfo)), ]


varPart <- fitExtractVarPartModel(exp2, form, bsinfo)

plotVarPart(varPart)

vp <- sortCols(varPart)

plotPercentBars(vp)

vp


#### remove medium
mdata <- merged_data[,1:12]
mdata <- mdata[mdata$stimulation !="M",]


merged_data2 <- merged_data[merged_data$stimulation !="M",]

mdatam <- merged_data2[,c("age_yrs", "gender","group", )]
mdatam <- cbind(mdata, mdatam)

mdatam <- na.omit(mdatam)

mod1 <- lmer( IL_8_1 ~ stimulation + age_yrs + gender + group + (1|ID),
              
              data = mdatam
)


library(car)
anova_results <- Anova(mod1, type = 2)  # Type II SS (recommended)
print(anova_results)


anova_results$`Sum Sq` / sum(anova_results$`Sum Sq`)

## the problem is random factor for variancepartion model, linear model can not taken the random
#into account, chatgpt said could use partR2()