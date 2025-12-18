rm(list = ls())

# Load required libraries
library(esvis)  # Load library for effect size visualization
library(pheatmap)  # Load library for creating heatmaps
library(tidyverse)  # Load library for data manipulation and visualization
library(dplyr)

#
Luminex_CAP_COVID_HV <- read.csv("Original_data/Luminex_CAP_COVID_HV_log.csv")

studydata_patients <- read.csv("Original_data/studydata_patients.csv")
studydata_HV <- read.csv("Original_data/studydata_HV.csv")
studydata <- rbind(studydata_patients, studydata_HV)

data <- merge(Luminex_CAP_COVID_HV, studydata, 
              by.x = "ID", by.y = "EB_id")

data_1 <- data[,c("CCL2","CCL4","IL_1RA","IL_8_1","TNF","CCL3","IL_1beta","IL_6","IL_10","stimulation","group")]
markers <- c("CCL2","CCL4","IL_1RA","IL_8_1","TNF","CCL3","IL_1beta","IL_6","IL_10")

data_long <- data_1 %>%
  pivot_longer(
    cols      = all_of(markers),
    names_to  = "marker",
    values_to = "value"
  )
data_long$marker <- factor(
  data_long$marker,
  levels = c("CCL2", "CCL3", "CCL4", "TNF", "IL_6", "IL_8_1", "IL_1beta", "IL_10", "IL_1RA")
)
data_long$stimulation <- factor(
  data_long$stimulation,
  levels = c("M", "LPS", "Kpneu", "Spneu")
)

#Reorder factor levels so HV is first, CAP second
data_long$group <- factor(data_long$group, levels = c("HV", "CAP"))

#Plot with manual color assignment
ggplot(data_long, aes(x = stimulation, y = value, fill = group)) +
  geom_boxplot() +
  facet_wrap(~ marker, ncol = 3) +
  scale_fill_brewer(palette = "Set2") +
  theme_bw(base_size = 16) +
  # Make facet labels (the text on top) smaller
  theme(
    strip.text = element_text(size = 10),      # facet (panel) label size
    plot.title  = element_text(size = 14)      # main title size
  ) +
  labs(
    x = "Stimulation",
    y = "Marker Value (log10 transformed)",
    title = "Profile of whole blood stimulation raw data"
  )



#minis####
data_2 <- data[,c("minis_M_CCL2","minis_M_CCL4","minis_M_IL_1RA","minis_M_IL_8_1","minis_M_TNF","minis_M_CCL3","minis_M_IL_1beta","minis_M_IL_6","minis_M_IL_10","stimulation","group")]
data_2 <- data_2 %>% filter(!stimulation == "M")
markers <- c("minis_M_CCL2","minis_M_CCL4","minis_M_IL_1RA","minis_M_IL_8_1","minis_M_TNF","minis_M_CCL3","minis_M_IL_1beta","minis_M_IL_6","minis_M_IL_10")

data_long <- data_2 %>%
  pivot_longer(
    cols      = all_of(markers),
    names_to  = "marker",
    values_to = "value"
  )

data_long$marker <- factor(
  data_long$marker,
  levels = c("minis_M_CCL2", "minis_M_CCL3", "minis_M_CCL4", "minis_M_TNF", "minis_M_IL_6", "minis_M_IL_8_1", "minis_M_IL_1beta", "minis_M_IL_10", "minis_M_IL_1RA")
)
data_long$stimulation <- factor(
  data_long$stimulation,
  levels = c("M", "LPS", "Kpneu", "Spneu")
)

#Reorder factor levels so HV is first, CAP second
data_long$group <- factor(data_long$group, levels = c("HV", "CAP"))

#Plot with manual color assignment
ggplot(data_long, aes(x = stimulation, y = value, fill = group)) +
  geom_boxplot() +
  facet_wrap(~ marker, ncol = 3) +
  scale_fill_brewer(palette = "Set2") +
  theme_bw(base_size = 16) +
  # Make facet labels (the text on top) smaller
  theme(
    strip.text = element_text(size = 10),      # facet (panel) label size
    plot.title  = element_text(size = 14)      # main title size
  ) +
  labs(
    x = "Stimulation",
    y = "Marker Value (log10 transformed)",
    title = "Profile of whole blood stimulation minus medium"
  )
