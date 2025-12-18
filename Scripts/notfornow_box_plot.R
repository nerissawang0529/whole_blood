rm(list = ls())

#Chemokines: CCL2, CCL3, CCL4
#Pro-Inflammatory Cytokines: IL-1β, IL-6, TNF, IL-8
#Anti-Inflammatory Cytokines: IL-1RA, IL-10

#ready data
studydata <- read.csv("Original_data/studydata.csv")
luminex_updated_3_with_m <- read.csv("Original_data/luminex_updated_3_with_m.csv")

#combine CAP,COVID-19,HV group
luminex_studydata <- merge(luminex_updated_3_with_m, 
                           studydata[, c("EB_id", "group")], 
                        by.x = "ID", by.y = "EB_id", 
                        all.x = TRUE)

luminex_studydata_without_M <- filter(luminex_studydata, stimulation != 'M')

#box plot
# Load necessary libraries
library(ggplot2)
library(dplyr)


# Create a combined factor for group and stimulation for plotting
luminex_studydata_without_M <- luminex_studydata_without_M %>%
  mutate(combined_group = interaction(group, stimulation, sep = "_"))

# Plot the data
library(ggplot2)

library(ggplot2)

p <- ggplot(luminex_studydata_without_M, aes(x = combined_group, y = minis_M_IL_10, fill = group)) +
  geom_boxplot(outlier.shape = NA) + # Suppress the black outliers from the boxplot
  geom_jitter(aes(fill = group), shape = 21, color = "black", position = position_jitter(width = 0.2, height = 0), alpha = 0.7, stroke = 0.5) + # Add jitter with black frame and group color fill
  scale_x_discrete(
    labels = c("HV_PP", "HV_KP", "HV_LPS", "CAP_PP", "CAP_KP", "CAP_LPS", "COVID_PP", "COVID_KP", "COVID_LPS"),
    limits = c("HV_PP", "HV_KP", "HV_LPS", "CAP_PP", "CAP_KP", "CAP_LPS", "COVID_PP", "COVID_KP", "COVID_LPS")
  ) +
  scale_fill_manual(values = c("HV" = "#a5d6a7", "CAP" = "#fbc02d", "COVID" = "#ce93d8")) + # Sage green for HV, mustard yellow for CAP, light lavender for COVID
  coord_cartesian(ylim = c(-800, 5000)) +  
  labs(x = "Group and Stimulation", y = "IL-10 minis medium (pg/ml)") +
  theme_minimal(base_size = 15) + # Increase base font size for larger text
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12), # Larger x-axis text
    axis.text.y = element_text(size = 12), # Larger y-axis text
    axis.title = element_text(size = 14), # Larger axis titles
    legend.title = element_text(size = 14), # Larger legend title
    legend.text = element_text(size = 12) # Larger legend text
  )

p


luminex_studydata_without_M$minis_M_IL_10


##export figure
file_path <- "Documents/Luminex/R_code/Figure/IL-10.svg" 
ggsave(file_path, plot = p, width = 8, height = 8)

#CCL2 -100000, 150000; CCL4 -50000, 300000; IL-1RA -150000, 300000; IL-8 -250000, 700000;
#TNF -10000, 50000;    CCL3 -80000, 350000; IL-1β  -20000, 100000;  IL-6 -20000, 130000;
#IL-10 -800, 5000

##export form
destination_folder <- "Documents/Luminex/R_code/Original_data/" 
export_file_name <- "luminex_studydata_without_M.csv" 
write.csv(luminex_studydata_without_M, file.path(destination_folder, export_file_name), row.names = FALSE)
cat("File exported to:", file.path(destination_folder, export_file_name), "\n")
#normal distrubution or non-normal distrubution