#rm(list = ls())

#Chemokines: CCL2, CCL3, CCL4
#Pro-Inflammatory Cytokines: IL-1β, IL-6, TNF, IL-8
#Anti-Inflammatory Cytokines: IL-1RA, IL-10

#ready data
studydata_patients <- read.csv("Original_data/studydata_patients.csv")
#studydata_HV <- read.csv("Original_data/studydata_HV.csv")
#studydata_combined <- rbind(studydata_patients, studydata_HV)

Luminex_CAP_COVID_HV <- read.csv("Original_data/Luminex_CAP_COVID_HV.csv")
merged_data <- merge(Luminex_CAP_COVID_HV, studydata_patients, 
                     by.x = "ID", by.y = "EB_id")

#box plot
# Load necessary libraries
library(ggplot2)
library(dplyr)


# Load the necessary library
library(ggplot2)

# Assuming 'luminex_studydata_without_M' is your dataframe
# Create the box plot
# Reorder and remove NA in the pathogen column
merged_data$pathogen <- factor(merged_data$pathogen, 
                                    levels = c("aureus", "haemo", "influenza", "Strep", "other", "mixed_infection", "No_pathogen", "Not"))

# Plotting code
p <- ggplot(merged_data, aes(x = pathogen, y = minis_M_IL_6, fill = stimulation)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75)) + # Position boxes with dodge
  geom_jitter(shape = 21, color = "black", position = position_dodge(width = 0.75), alpha = 0.3, stroke = 0.3) + # Use dodge for jitter as well
  labs(x = "Pathogen", y = "IL-6 Concentration") +
  scale_fill_manual(values = c("PP" = "#ADD8E6", "KP" = "#FFDAB9", "LPS" = "#BC8F8F")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(fill = "Stimulation") +
  coord_cartesian(ylim = c(-20000, 130000)) +  
  labs(x = "Stimulation", y = "IL-6 minis medium (pg/ml)") +
  theme_minimal(base_size = 15) + # Increase base font size for larger text
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12), # Larger x-axis text
    axis.text.y = element_text(size = 12), # Larger y-axis text
    axis.title = element_text(size = 14), # Larger axis titles
    legend.title = element_text(size = 14), # Larger legend title
    legend.text = element_text(size = 12) # Larger legend text
  )

# Display plot
print(p)


#to get the different pathogen numbers in high and low il-6 group
#need to run the what is wrong related code
library(dplyr)

# 
hil6_data <- luminex_long %>%
  dplyr::select(ID, HiIL6) %>%
  distinct()
merged_data <- merged_data %>%
  mutate(ID = as.character(ID)) %>%
  left_join(hil6_data, by = "ID")

merged_data <- merged_data %>%
  left_join(hil6_data, by = "ID")
merged_unique <- merged_data %>%
  distinct(ID, .keep_all = TRUE)


merged_unique %>%
  dplyr::count(pathogen, HiIL6.x) %>%
  group_by(pathogen) %>%
  mutate(prop = n / sum(n)) %>%
  pivot_wider(names_from = HiIL6.x, values_from = prop, values_fill = 0)
library(dplyr)
library(tidyr)

