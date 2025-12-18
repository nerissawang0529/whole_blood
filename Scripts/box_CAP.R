rm(list = ls())

#Chemokines: CCL2, CCL3, CCL4
#Pro-Inflammatory Cytokines: IL-1β, IL-6, TNF, IL-8
#Anti-Inflammatory Cytokines: IL-1RA, IL-10

#ready data
studydata_patients <- read.csv("Original_data/studydata_patients.csv")

Luminex_CAP_COVID_HV <- read.csv("Original_data/Luminex_CAP_COVID_HV_log.csv")
merged_data <- merge(Luminex_CAP_COVID_HV, studydata_patients, 
                     by.x = "ID", by.y = "EB_id")

num_markers <- c("CCL2","CCL4","IL_1RA","IL_8_1","TNF","CCL3","IL_1beta","IL_6","IL_10")

merged_data[num_markers] <- lapply(merged_data[num_markers], function(x) as.numeric(as.character(x)))

# **Scale the numeric cytokine markers**
# or log the data
#merged_data[num_markers] <- scale(merged_data[num_markers])
#merged_data[, num_markers] <- log10(merged_data[, num_markers])
#box plot
# Load necessary libraries
library(ggplot2)
library(dplyr)


# Assuming 'luminex_studydata_without_M' is your dataframe
# Create the box plot
# Reorder and remove NA in the pathogen column
#luminex_combined$pathogen <- factor(luminex_combined$pathogen, 
#                                    levels = c("aureus", "haemo", "influenza", "Strep", "other", "mixed_infection", "No_pathogen", "Not"))

# Plotting code
data_long <- merged_data %>%
  pivot_longer(cols = all_of(num_markers), 
               names_to = "marker", 
               values_to = "value")
data_long$marker <- factor(
  data_long$marker,
  levels = c("CCL2", "CCL3", "CCL4", "TNF", "IL_6", "IL_8_1", "IL_1beta", "IL_10", "IL_1RA")
)
data_long$stimulation <- factor(
  data_long$stimulation,
  levels = c("M", "LPS", "Kpneu", "Spneu")
)
ggplot(data_long, aes(x = stimulation, y = value, fill = stimulation)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75)) + 
  geom_jitter(shape = 21, color = "black", position = position_dodge(width = 0.75), 
              alpha = 0.3, stroke = 0.3) + 
  facet_wrap(~ marker, ncol = 3, scales = "free_y") +  
  scale_fill_brewer(palette = "Set2") +  
  theme_bw(base_size = 16) +  
  theme(
    strip.text = element_text(size = 10),  
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12), 
    axis.text.y = element_text(size = 12), 
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 12)
  ) +
  labs(
    x = "Stimulation",
    y = "Marker Value  (log10 transformed)",
    title = "Profile of CAP patients"
  )










studydata <- read.csv("Original_data/studydata_patients.csv")
studydata$EB_id

studydata_all <- read.csv("Original_data/studydata.csv")
luminex_updated_3_with_m <- read.csv("Original_data/luminex_updated_3_with_m.csv")


#combine CAP,HV group
luminex_studydata <- merge(luminex_updated_3_with_m, 
                           studydata[, c("EB_id", "group","pathogen")], 
                           by.x = "ID", by.y = "EB_id", 
                           all.x = TRUE)
#export clinical_marker_unique data
destination_folder <- "Original_data/" 
export_file_name <- "Whole_blood_stimulation.csv" 
write.csv(luminex_studydata, file.path(destination_folder, export_file_name), row.names = FALSE)
cat("File exported to:", file.path(destination_folder, export_file_name), "\n")

luminex_studydata_without_M <- filter(luminex_studydata, stimulation != 'M')
luminex_studydata_without_M <- luminex_studydata_without_M %>% filter(group %in% c("CAP"))

#add HV
studydata_hv <- studydata_all %>% filter(group %in% c("HV"))
luminex_studydata_hv <- merge(luminex_updated_3_with_m, 
                              studydata_hv[, c("EB_id", "group")], 
                              by.x = "ID", by.y = "EB_id", 
                              all.x = TRUE)
luminex_studydata_hv_without_M <- filter(luminex_studydata_hv, stimulation != 'M')
luminex_studydata_hv_without_M <- luminex_studydata_hv_without_M %>% filter(group %in% c("HV"))
luminex_studydata_hv_without_M$pathogen <- "Not"


luminex_combined <- rbind(luminex_studydata_without_M, luminex_studydata_hv_without_M)
