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
luminex <- import("Original_data/Final_data_20250301.xlsx")


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


install.packages("stringr")
library(stringr)

# Columns to check for ">" values
biomarker_cols <- c("CCL2","CCL4","IL_1RA","IL_8_1","TNF","CCL3","IL_1beta","IL_6","IL_10")  # Add all biomarker columns where ">" may appear

# Function to update biomarker values
update_luminex <- function(df, biomarker_cols) {
  df <- df %>%
    arrange(ID, stimulation)  # Sort to ensure correct order of duplicates
  
  for (col in biomarker_cols) {
    # Find rows with the same ID and stimulation, where the biomarker has a ">" symbol
    duplicates <- df %>%
      filter(str_detect(!!sym(col), ">")) %>%
      dplyr::select(ID, stimulation) %>%
      unique()
    
    for (i in 1:nrow(duplicates)) {
      # Filter rows for this specific ID and stimulation
      rows <- which(df$ID == duplicates$ID[i] & df$stimulation == duplicates$stimulation[i])
      
      # Check if there are two rows (the duplicate case)
      if (length(rows) == 2) {
        # Replace ">" value in the first row with the numeric value from the second row
        replacement_value <- df[[col]][rows[2]]
        
        # Update the first row's value
        df[[col]][rows[1]] <- replacement_value
      }
    }
  }
  
  return(df)
}

# Apply the function to the luminex dataframe
library(dplyr)
luminex_updated <- update_luminex(luminex, biomarker_cols)

# Remove duplicates, keeping only the first row for each ID and stimulation
luminex_updated_2 <- luminex_updated %>%
  group_by(ID, stimulation) %>%
  slice(1) %>%
  ungroup()



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
  for (i in 2:10) {  # Iterate over biomarker columns
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

#> and space

luminex_updated_3 <- as.data.frame(lapply(new_data, function(x) {
  sapply(x, function(cell) {
    if (grepl(">", cell) || cell == "" || is.na(cell) || tolower(cell) == "blank") {  
      NA  # Replace with NA if '>', empty, or 'blank' is found
    } else {
      cell  # Leave cell as is if none of these conditions are met
    }
  })
}))

str(luminex_updated_3)

#into number 
luminex_updated_3[biomarker_cols] <- lapply(luminex_updated_3[biomarker_cols], function(x) as.numeric(as.character(x)))
luminex_updated_3$stimulation[luminex_updated_3$stimulation == "PP"] <- "Spneu"
luminex_updated_3$stimulation[luminex_updated_3$stimulation == "KP"] <- "Kpneu"

#into log then export as log data####
# Apply log10 transformation
luminex_updated_3 <- luminex_updated_3 %>%
  mutate(across(all_of(biomarker_cols), ~ log10(.), .names = "{.col}"))

#minis M
# Load necessary libraries

# Step 1: Create a dataframe of 'M' values for each biomarker and ID
m_values <- luminex_updated_3 %>%
  filter(stimulation == "M") %>%
  dplyr::select(ID, all_of(biomarker_cols)) %>%
  rename_with(~ paste0("M_", .), -ID)

# Step 2: Join the 'M' values back to the original dataframe by ID
luminex_updated_3_with_m <- luminex_updated_3 %>%
  left_join(m_values, by = "ID")

# Step 3: For each biomarker, subtract the corresponding 'M' value if stimulation is KP, LPS, or PP
luminex_updated_3_with_m <- luminex_updated_3_with_m %>%
  mutate(across(all_of(biomarker_cols), 
                ~ if_else(stimulation %in% c("Kpneu", "LPS", "Spneu"), . - get(paste0("M_", cur_column())), NA_real_), 
                .names = "minis_M_{col}"))

# Step 4: Remove temporary columns if needed
luminex_updated_3_with_m <- luminex_updated_3_with_m %>%
  dplyr::select(-starts_with("M_"))

#into number 
luminex_updated_3_with_m[biomarker_cols] <- lapply(luminex_updated_3_with_m[biomarker_cols], function(x) as.numeric(as.character(x)))



#i used the luminex_updated_3_with_m$IL_8_1 is NA for not retested IDs for now.
#not_retested_COVID <- luminex_updated_3_with_m$ID[is.na(luminex_updated_3_with_m$IL_8_1)]
#the ID is 3024, 1194, 1195, 1210, 1227, 1228, 3095, 3098, 3099, 3104, 3115, 3128, 3165, 3166, 3167, 3181, 3186, 3201, 3207, 3208, 3218, 3219, 3227, 3245, 3286, 1149, 1160, 1186, 3174, 3291, 3248
ids_to_remove <- c(3024, 1194, 1195, 1210, 1227, 1228, 3095, 3098, 3099, 3104, 
                   3115, 3128, 3165, 3166, 3167, 3181, 3186, 3201, 3207, 3208, 
                   3218, 3219, 3227, 3245, 3286, 1149, 1160, 1186, 3174, 3291, 3248)

luminex_updated_3_with_m <- luminex_updated_3_with_m %>%
  filter(!ID %in% ids_to_remove)

#exclude which are not penumounia 
notCOVID_notCAP <- read.csv("Original_data/notCOVID_notCAP.csv")
luminex_updated_3_with_m <- luminex_updated_3_with_m %>%
  filter(!ID %in% notCOVID_notCAP$EB_id)

#exclude which don't have M value (too high and not retested)
ids_to_remove_2 <- c(1069,1004)
luminex_updated_3_with_m <- luminex_updated_3_with_m %>%
  filter(!ID %in% ids_to_remove_2)

Luminex_CAP_COVID_HV <- luminex_updated_3_with_m



##export form
destination_folder <- "Original_data/" 
export_file_name <- "Luminex_CAP_COVID_HV_log.csv"  #log markers then substract medium
write.csv(Luminex_CAP_COVID_HV, file.path(destination_folder, export_file_name), row.names = FALSE)
cat("File exported to:", file.path(destination_folder, export_file_name), "\n")
