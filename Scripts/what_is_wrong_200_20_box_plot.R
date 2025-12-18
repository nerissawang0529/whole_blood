#rm(list = ls())

#ready for studydata data####

#IDs should be noticed
#Stijn: for patient 1090 the tubes where somehow unlabeled, keep it out for now
#Erik:  can use:    2121, 2125, 2116  
#       cannot use: 2120, 3345, 3137 # 2120 and 3345 are not in the excel form
## install.packages("pacman")
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rxio, ggplot2, naniar,
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
library(readxl)

luminex <- read_excel("Original_data/Final_data_20250320.xlsx")
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

library(dplyr)
library(ggplot2)

#show all markers
# Load required libraries
library(dplyr)
library(ggplot2)
library(purrr)

# Define markers of interest
markers_to_process <- c("CCL2", "CCL3", "CCL4", "IL_6", "IL_1beta", "IL_10", "IL_1RA", "TNF")

# Function to process a single marker
process_marker <- function(marker_name) {
  message("Processing: ", marker_name)
  
  tnf_raw <- luminex_long %>%
    filter(marker == marker_name, Dilution %in% c(2, 200)) %>%
    filter(!grepl(">", value)) %>%
    mutate(value = as.numeric(value))
  
  if (nrow(tnf_raw) == 0) {
    warning("No valid rows for marker: ", marker_name)
    return(tibble())  # return empty tibble so map_dfr still works
  }
  
  valid_ids <- tnf_raw %>%
    group_by(ID_stimulation) %>%
    summarise(n_dilutions = n_distinct(Dilution), .groups = "drop") %>%
    filter(n_dilutions == 2) %>%
    pull(ID_stimulation)
  
  filtered <- tnf_raw %>%
    filter(ID_stimulation %in% valid_ids) %>%
    mutate(
      Dilution_group = case_when(
        Dilution == 2 ~ "2x",
        Dilution == 200 ~ "200x"
      )
    )
  
  estimated_20x <- filtered %>%
    filter(Dilution == 200) %>%
    mutate(
      value = value / 10,
      Dilution_group = "200x_as_20x"
    )
  
  result <- bind_rows(filtered, estimated_20x) %>%
    mutate(marker = marker_name)
  
  return(result)
}


# Process all markers and combine
all_markers_data <- map_dfr(markers_to_process, process_marker)

# Set factor order for plotting
all_markers_data$Dilution_group <- factor(all_markers_data$Dilution_group,
                                          levels = c("2x", "200x", "200x_as_20x"))

# Loop through each marker and display plot
for (m in markers_to_process) {
  p <- all_markers_data %>%
    filter(marker == m) %>%
    ggplot(aes(x = Dilution_group, y = value, fill = Dilution_group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.15, size = 0.6, alpha = 0.3) +
    scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000)) +
    facet_wrap(~ group) +
    labs(
      title = paste("Cytokine:", m),
      x = "Dilution group",
      y = "Concentration (log10)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      strip.text = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(angle = 0),
      panel.grid.minor = element_blank()
    )
  
  # Show the plot in RStudio Viewer
  print(p)
}
