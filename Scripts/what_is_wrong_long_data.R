rm(list = ls())

#ready for studydata data####

#IDs should be noticed
#Stijn: for patient 1090 the tubes where somehow unlabeled, keep it out for now
#Erik:  can use:    2121, 2125, 2116  
#       cannot use: 2120, 3345, 3137 # 2120 and 3345 are not in the excel form
## install.packages("pacman")
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
               tableone, kmed, survival, survminer, skimr, lattice, reshape2, grDevices, ggrepel)

install.packages("dplyr")
library(dplyr)
library(MASS)
library(magrittr)
library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)
library(readxl)

#####get ready for the luminex data####
luminex <- import("Original_data/Final_data_20250320.xlsx")
luminex <- luminex %>%
  mutate(Kit_name = gsub("plate", "", Kit_name))
studydata_patients <- read.csv("Original_data/studydata_patients.csv")
studydata_HV <- read.csv("Original_data/studydata_HV.csv")
library(dplyr)
library(dplyr)


#to get the id, i added _N for the 1016 plate3456 after the ids in the 'Final_data_20241026'.
luminex$Modified_Position <- ifelse(
  grepl("_[A-Z]+_[A-Z]$", luminex$Position), 
  sub(".*_(\\d+_[A-Z]+)_.*", "\\1", luminex$Position),            # Case like '3328_KP' or '1026_LPS'
  ifelse(grepl("_[A-Z]$", luminex$Position), 
         sub(".*_(\\d+)_.*", "\\1_M", luminex$Position),           # Case like '1173_Y' or '2036_N' -> '1173_M' or '2036_M'
         ifelse(grepl("_\\d+$", luminex$Position),
                sub(".*_(\\d+)$", "\\1_M", luminex$Position),      # Case like '2036' (ends in a number) -> '2036_M'
                luminex$Position))                                 # Leave as-is if no match
)
luminex$Modified_Position_2 <- gsub(".*\\s", "", luminex$Modified_Position)
luminex$Modified_Position_3 <- gsub("_N", "", luminex$Modified_Position_2)
luminex$ID_stimulation <-luminex$Modified_Position_3
luminex <- luminex %>%
  separate(Modified_Position_3, into = c("ID", "stimulation"), sep = "_")
luminex$Modified_Position <- NULL
luminex$Modified_Position_2 <- NULL
luminex$Position <- NULL
#one sitimulaiton is wrong, delete this sample
luminex <- luminex[luminex$ID != 2098, ]

#merge the group into dataframe
group_data_HV <- studydata_HV %>%
  dplyr::select(EB_id, group) 

group_data_CAP <- studydata_patients %>%
  dplyr::select(EB_id, group) 

group_data <- bind_rows(group_data_HV, group_data_CAP)
names(group_data)[names(group_data) == "EB_id"] <- "ID"

luminex <- luminex %>%
  mutate(ID = as.character(ID))

group_data <- group_data %>%
  mutate(ID = as.character(ID))

luminex <- luminex %>%
  left_join(group_data, by = "ID") %>%
  mutate(group = ifelse(is.na(group), "COVID", group))  # Assign "COVID" to missing groups


#the ID is 3024, 1194, 1195, 1210, 1227, 1228, 3095, 3098, 3099, 3104, 3115, 3128, 3165, 3166, 3167, 3181, 3186, 3201, 3207, 3208, 3218, 3219, 3227, 3245, 3286, 1149, 1160, 1186, 3174, 3291, 3248
ids_to_remove <- c(3024, 1194, 1195, 1210, 1227, 1228, 3095, 3098, 3099, 3104, 
                   3115, 3128, 3165, 3166, 3167, 3181, 3186, 3201, 3207, 3208, 
                   3218, 3219, 3227, 3245, 3286, 1149, 1160, 1186, 3174, 3291, 3248)

luminex <- luminex %>%
  filter(!ID %in% ids_to_remove)

#exclude which are not penumounia 
notCOVID_notCAP <- read.csv("Original_data/notCOVID_notCAP.csv")
luminex <- luminex %>%
  filter(!ID %in% notCOVID_notCAP$EB_id)

#exclude which don't have M value (too high and not retested)
ids_to_remove_2 <- c(1069,1004)
luminex <- luminex %>%
  filter(!ID %in% ids_to_remove_2)

luminex_updated_2 <- luminex

# replace the '<'
luminex_updated_2 <- as.data.frame(lapply(luminex_updated_2, function(x) {
  sapply(x, function(cell) {
    if (grepl("\\*", cell)) {  # Check if * is in the cell
      gsub("\\*", "", cell)    # Remove * symbol if it appears
    } else {
      cell  # Leave cell as is if no * symbol is found
    }
  })
}))
luminex_updated_2 <- as.data.frame(lapply(luminex_updated_2, function(x) {
  sapply(x, function(cell) {
    if (grepl("\\,", cell)) {  
      gsub(",", ".", cell)    
    } else {
      cell  # Leave cell as is if no * symbol is found
    }
  })
}))


tempt_data_list <- luminex_updated_2 %>%
  group_split(Day) 

find_non_numeric_elements <- function(vec) {
  # Attempt to convert each element and check for NAs
  non_numeric_indices <- which(is.na(suppressWarnings(as.numeric(vec))))
  
  # Return the non-numeric elements
  return(vec[non_numeric_indices])
}


new_data_list <- NULL
for (j in 1:length(tempt_data_list)) {  
  a <- tempt_data_list[[j]]  # Select data for Day j
  for (i in 4:12) {  # Iterate over biomarker columns
    cur_column <- a[, i] %>% unlist()
    # Find the minimum numeric value within the SAME Day
    min_value <- cur_column[which(!is.na(suppressWarnings(as.numeric(cur_column))))] %>% 
      as.numeric %>% 
      min(., na.rm = TRUE)
    # Replace "<" values within the same Day using the minimum for that Day
    cur_column[grepl("<", cur_column)] <- min_value / sqrt(2)
    # Update the column in the current day's dataframe
    a[, i] <- cur_column
  }
  # Store the modified data frame in the new list
  new_data_list[[j]] <- a
}
# Combine all data frames in the list into a single data frame
new_data <- do.call(rbind, new_data_list)
luminex <- new_data



#make the data into long
library(dplyr)
library(tidyr)

# Define columns to keep
keep_cols <- c("Day", "Kit_name", "Dilution", "ID", "stimulation", "ID_stimulation", "group")

# Convert to long format
luminex_long <- luminex %>%
  pivot_longer(cols = -all_of(keep_cols),  # Pivot all columns except the kept ones
               names_to = "marker",        # New column for marker names
               values_to = "value")        # New column for marker values

# >#####
library(dplyr)

# Convert value to character (if it's not already)
luminex_long <- luminex_long %>%
  mutate(value = as.character(value))

# Identify rows where value contains ">"
greater_than_rows <- luminex_long %>%
  filter(grepl(">", value))  # Finds rows with ">"

# Identify rows WITHOUT ">" to check for duplicates
non_greater_than_rows <- luminex_long %>%
  filter(!grepl(">", value))

# Find IDs that have both ">" and normal values
common_ids <- greater_than_rows %>%
  dplyr::select(ID, stimulation, marker) %>%
  distinct()

# Replace ">" rows with their corresponding non ">" rows
filtered_non_greater_than_rows <- non_greater_than_rows %>%
  semi_join(common_ids, by = c("ID", "stimulation", "marker"))

# Keep rows where ">" is not present AND "Dilution != 200"
remaining_rows <- non_greater_than_rows %>%
  anti_join(common_ids, by = c("ID", "stimulation", "marker")) %>%
  filter(Dilution != 200)

# Combine the replacements with valid remaining rows
luminex_long <- bind_rows(filtered_non_greater_than_rows, remaining_rows)

# Ensure uniqueness
luminex_long <- luminex_long %>%
  distinct()


#check for HV
# filter data in one step
luminex_M_CCL2_HV <- luminex_long %>%
  filter(stimulation == "M", marker == "CCL2", group == "HV") %>%
  mutate(
    value = as.numeric(value),
    Dilution = as.factor(Dilution),
    Day = as.factor(Day),
    Kit_name = as.factor(Kit_name),
    level = log10(value)  # Log-transform value
  )
hist(luminex_M_CCL2_HV$level)

# fit Linear Model
model <- lm(level ~ Dilution + Kit_name, data = luminex_M_CCL2_HV)

# rint results
summary(model)
anova(model)

# boxplot
ggplot(luminex_M_CCL2_HV, aes(x = Kit_name, y = level, fill = as.factor(Dilution))) +
  geom_boxplot() +
  labs(title = "CCL2 Levels by Kit (non-infection)", 
       y = "Log10 CCL2", 
       x = "Kit Name", 
       fill = "Dilution") +  # Legend label
  theme_minimal()

#check for CAP
# filter data in one step
luminex_M_IL_1beta_CAP <- luminex_long %>%
  filter(stimulation == "M", marker == "IL_1beta", group == "CAP") %>%
  mutate(
    value = as.numeric(value),
    Dilution = as.factor(Dilution),
    Day = as.factor(Day),
    Kit_name = as.factor(Kit_name),
    level = log10(value)  # Log-transform value
  )
hist(luminex_M_IL_1beta_CAP$level)
# fit Linear Model
model <- lm(level ~ Dilution + Kit_name, data = luminex_M_IL_1beta_CAP)

# print results
summary(model)
anova(model)

# goxplot
ggplot(luminex_M_IL_1beta_CAP, aes(x = Kit_name, y = level, fill = as.factor(Dilution))) +
  geom_boxplot() +
  labs(title = "CCL2 Levels by Kit (CAP)", 
       y = "Log10 CCL2", 
       x = "Kit Name", 
       fill = "Dilution") +  # Legend label
  theme_minimal()



#check for COVID
luminex_M_CCL2_COVID <- luminex_long %>%
  filter(stimulation == "M", marker == "CCL2", group == "COVID") %>%
  mutate(
    value = as.numeric(value),
    Dilution = as.factor(Dilution),
    Day = as.factor(Day),
    Kit_name = as.factor(Kit_name),
    level = log10(value)  # Log-transform value
  )
hist(luminex_M_CCL2_COVID$level)

# fit Linear Model
model <- lm(level ~ Dilution + Kit_name, data = luminex_M_CCL2_COVID)

# print results
summary(model)
anova(model)

# goxplot
ggplot(luminex_M_CCL2_COVID, aes(x = Kit_name, y = level, fill = as.factor(Dilution))) +
  geom_boxplot() +
  labs(title = "CCL2 Levels by Kit (COVID)", 
       y = "Log10 CCL2", 
       x = "Kit Name", 
       fill = "Dilution") +  # Legend label
  theme_minimal()

