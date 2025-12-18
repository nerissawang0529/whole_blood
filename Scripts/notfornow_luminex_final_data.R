
#ready for studydata data####

#IDs should be noticed
#Stijn: for patient 1090 the tubes where somehow unlabeled, keep it out for now
#Erik:  can use:    2121, 2125, 2116  
#       cannot use: 2120, 3345, 3137 # 2120 and 3345 are not in the excel form

## install.packages("pacman")
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
               tableone, kmed, survival, survminer, skimr, lattice, reshape2, grDevices, ggrepel)


#ready the clinical data
studydata_luminex <- import("Documents/Luminex/R_code/Original_data/studydata.csv")

#####get ready for the luminex data####
luminex <- import("Original_data/Final_data_20241026.xlsx")
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




#used the data only > for 200 times dilution
# Load necessary libraries
library(dplyr)
library(stringr)

# Define biomarker columns where ">" may appear
biomarker_cols <- c("CCL2", "CCL4", "IL_1RA", "TNF", "CCL3", "IL_6", "IL_10", "IL_8_1", "IL_1beta")

# Function to update biomarker values
update_luminex <- function(df, biomarker_cols) {
  df <- df %>%
    arrange(ID, stimulation)  # Sort to ensure correct order of duplicates
  
  for (col in biomarker_cols) {
    # Find rows with the same ID and stimulation, where the biomarker has a ">" symbol
    duplicates <- df %>%
      filter(str_detect(.[[col]], ">")) %>%  # Directly access column by name
      distinct(ID, stimulation)  # Unique combinations of ID and stimulation
    
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
luminex_updated <- update_luminex(luminex, biomarker_cols)


# Remove duplicates, keeping only the first row for each ID and stimulation
luminex_updated_2 <- luminex_updated %>%
  group_by(ID, stimulation) %>%
  slice(1) %>%
  ungroup()








luminex_updated_3 <- as.data.frame(lapply(luminex_updated_2, function(x) {
  sapply(x, function(cell) {
    if (grepl("\\*", cell)) {  
      gsub("\\*", "", cell)    
    } else {
      cell  
    }
  })
}))
luminex_updated_4 <- as.data.frame(lapply(luminex_updated_3, function(x) {
  sapply(x, function(cell) {
    if (grepl("\\,", cell)) {  
      gsub(",", ".", cell)    
    } else {
      cell 
    }
  })
}))
luminex_updated_5 <- as.data.frame(lapply(luminex_updated_4, function(x) {
  sapply(x, function(cell) {
    if (grepl("\\>", cell)) {  
      gsub(">", "NA", cell)    
    } else {
      cell  
    }
  })
}))


##export form
destination_folder <- "Documents/Luminex/R_code/Original_data/" 
export_file_name <- "luminex_updated_2.csv" 
write.csv(luminex_updated_2, file.path(destination_folder, export_file_name), row.names = FALSE)
cat("File exported to:", file.path(destination_folder, export_file_name), "\n")

# deal with '<' by Gpt
#in the form, there is a column called 'Day', there is 0903_plate1_2 ,0903_plate3_4 ,0903_plate5_6 ,0904_plate1_2 ,0904_plate3_4, 0904_plate5_6, 1016_plate1_2, 1016_plate3_4, 1016_plate5_6 ,Day_1, Day_2,Day_3, kit4, kit5, kit7. 
#and there are some markers: "CCL2" ,"CCL4" ,"IL_1RA","IL_8_1" ,"TNF", "CCL3" ,"IL_1beta" ,"IL_6", "IL_10" . 
#Now please help me replace the '<', using the rules: within that Day, using 'divide the lowest value of that day by the square root of 2' of to replace the '<'. For example, when the 'Day' is 0903_plate1_2, in column "CCL2", there is a '<', then change the '<' into 'divide the lowest value of that day by the square root of 2'

luminex_updated_2_gpt <- import("Documents/Luminex/R_code/Original_data/luminex_updata_2_gpt.csv")
luminex_updated_2_gpt$V1 <- NULL








luminex_updated_2[, 2:10] <- lapply(luminex_updated_2[, 2:10], function(x) {
    if(is.character(x)) {
        x <- gsub("\\*", "", x) # Remove asterisks
        x <- gsub(",", ".", x)  # Replace commas with dots
      }
  as.numeric(x) # Convert to numeric if possible
  })



luminex_updated <- merge(luminex_updated_2, studydata_luminex[, c("EB_id", "group")], 
                 by.x = "ID", by.y = "EB_id", all.x = TRUE)


# Apply the cleaning function only to columns 3 through 11



# Display the cleaned dataframe
print(luminex_updated)


#minis M
# Load necessary libraries
library(dplyr)

# Specify the biomarker columns explicitly
biomarker_cols <- c("CCL2", "CCL4", "IL_1RA", "TNF", "CCL3", "IL_6", "IL_10", "IL_8_1", "IL_1beta")

# Step 1: Create a dataframe of 'M' values for each biomarker and ID
m_values <- luminex_updated %>%
  filter(stimulation == "M") %>%
  select(ID, all_of(biomarker_cols)) %>%
  rename_with(~ paste0("M_", .), -ID)

# Step 2: Join the 'M' values back to the original dataframe by ID
luminex_updated_with_m <- luminex_updated %>%
  left_join(m_values, by = "ID")

# Step 3: For each biomarker, subtract the corresponding 'M' value if stimulation is KP, LPS, or PP
luminex_updated_with_m <- luminex_updated_with_m %>%
  mutate(across(all_of(biomarker_cols), 
                ~ if_else(stimulation %in% c("KP", "LPS", "PP"), . - get(paste0("M_", cur_column())), NA_real_), 
                .names = "minis_M_{col}"))

# Step 4: Remove temporary columns if needed
luminex_updated_with_m <- luminex_updated_with_m %>%
  select(-starts_with("M_"))

# Now luminex_updated_with_m should have new columns for each biomarker with 'minis_M' calculated correctly


#i didn't set the > <####


##export form
destination_folder <- "Documents/Luminex/R_code/Original_data/" 
export_file_name <- "luminex_marker_data.csv" 
write.csv(luminex, file.path(destination_folder, export_file_name), row.names = FALSE)
cat("File exported to:", file.path(destination_folder, export_file_name), "\n")


###for hist M
# Assuming your data is in a dataframe called 'data'
# Filter the data where stimulation is 'M'
filtered_data <- luminex_updated_with_m[luminex_updated_with_m$stimulation == "M", ]

# Plot the histogram for CCL2 in the filtered data
hist(filtered_data$IL_8_1, 
     main = "Histogram of IL_8_1 (Stimulation = M)", 
     xlab = "IL_8_1", 
     ylab = "Frequency", 
     col = "lightblue", 
     border = "black")
