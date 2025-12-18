#need to run what_is_wrong_long_data first to get the database


## install.packages("pacman")
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
               tableone, kmed, survival, survminer, skimr, lattice, reshape2, grDevices, ggrepel)

Bas <- readxl::read_xlsx("Original_data/Bas_figure1_data.xlsx")

library(dplyr)
library(tidyr)
library(stringr)

#clean data####
# Standardize column names (remove hidden spaces or special characters)
colnames(Bas) <- make.names(colnames(Bas))

Bas_long <- Bas %>%
  # Keep only rows where Inhoud starts with "WB_T"
  filter(str_detect(Inhoud, "^WB_T")) %>%
  
  # Split the ID column into ID and Day
  mutate(
    ID_raw = ID,
    ID = ifelse(str_detect(ID_raw, "_"), str_extract(ID_raw, "^[^_]+"), ID_raw),
    Day = ifelse(str_detect(ID_raw, "_"), str_extract(ID_raw, "[^_]+$"), "D0")
  ) %>%
  
  # Extract stimulation type from Inhoud column
  mutate(
    stimulation = case_when(
      str_detect(Inhoud, "M_sup") ~ "M",
      str_detect(Inhoud, "LPS_sup") ~ "LPS",
      str_detect(Inhoud, "KP_sup") ~ "KP",
      str_detect(Inhoud, "SP_sup") ~ "SP",
      TRUE ~ NA_character_
    )
  ) %>%
  
  # Pivot wide format to long format, assuming markers start from the IFN_gamma column
  pivot_longer(
    cols = starts_with("IFN_"):TNF_alpha,  # Alternatively: cols = 4:ncol(Bas)
    names_to = "marker",
    values_to = "value"
  ) %>%
  
  # Convert value column to numeric, handle "OOR <" by setting it to NA
  mutate(
    value = as.numeric(str_replace(value, "OOR <", NA_character_)),
    level = log10(value)
  ) %>%
  
  # Select and reorder columns
  dplyr::select(Day, ID, stimulation, marker, value, level)

# Show first few rows of the transformed data
head(Bas_long)

Bas_long_D0 <- Bas_long %>%
  filter(Day == "D0")


# add highIL6 group
# Identify high IL-6 IDs under M stimulation
highIL6 <- Bas_long_D0 %>%
  filter(marker == "IL_6", stimulation == "M", log10(value) > 3) %>%
  pull(ID) %>%
  unique()

# Add HiIL6 column to Bas_long
Bas_long_D0 <- Bas_long_D0 %>%
  mutate(
    HiIL6 = ifelse(ID %in% highIL6, "HighIL6_M", "LowIL6_M")
  )

##get the same ID in luminex_long
Bas_long_D0 <- Bas_long_D0 %>%
  filter(ID %in% unique(luminex_long$ID))

# Create unique ID-to-group mapping from luminex_long
id_group_lookup <- luminex_long %>%
  dplyr::select(ID, group) %>%
  distinct()

# Ensure ID is character in both dataframes (if needed)
Bas_long_D0$ID <- as.character(Bas_long_D0$ID)
id_group_lookup$ID <- as.character(id_group_lookup$ID)

# Safe merge based on ID only
Bas_long_D0 <- Bas_long_D0 %>%
  left_join(id_group_lookup, by = "ID")

#overview of the data of BAS
library(dplyr)
library(ggplot2)
library(ggbeeswarm)  # for geom_quasirandom

# Step 1: Define the marker list of interest
selected_markers <- c("CCL2", "CCL3", "CCL4", 
                      "IL_1β", "IL_6", "TNF_alpha", "IL_8", 
                      "IL_1RA", "IL_10")

# Step 2: Filter Bas_long_D0 to include only these markers
Bas_long_D0_filtered <- Bas_long_D0 %>%
  filter(marker %in% selected_markers)

# Step 3: Set stimulation order
Bas_long_D0_filtered <- Bas_long_D0_filtered %>%
  mutate(stimulation = factor(stimulation, levels = c("M", "LPS", "SP", "KP")))

# Step 4: Plot
plot_bas <- ggplot(data = Bas_long_D0_filtered, aes(x = stimulation, y = level)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_quasirandom(size = 0.8, aes(col = HiIL6, alpha = HiIL6)) +
  facet_grid(marker ~ group, scales = "free") +
  scale_color_manual(values = c("red", "green3")) +
  scale_alpha_manual(values = c(0.5, 0.8)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 10)
  ) +
  ggtitle("Xanthe & Bas: Individuals coloured red if IL6 > 3 in the medium condition") +
  ylab("log10(concentration [pg/ml])")

# Step 5: Show plot
plot_bas

ggsave(file="/Users/huiwang/Documents/Luminex/luminex_R/Figure/bas_il6.svg", plot=plot_bas, width=15, height=30)


#overview of the data Erik & Nerissa, same IDs as Bas
library(dplyr)
library(ggplot2)
library(ggbeeswarm)

# Subset luminex_long to include only IDs in Bas_long_D0_filtered
luminex_long$value <- as.numeric(luminex_long$value)

luminex_subset <- luminex_long %>%
  filter(ID %in% unique(Bas_long_D0_filtered$ID))

# Step 1: Define the selected markers
selected_markers <- c("IL_10", "IL_1RA", "IL_6", "TNF")

# Step 2: Filter data to keep only selected markers
luminex_subset_selected <- luminex_subset %>%
  filter(marker %in% selected_markers)

# Step 3: Set factor levels for stimulation to control x-axis order
luminex_subset_selected$stimulation <- factor(luminex_subset_selected$stimulation,
                                              levels = c("M", "LPS", "SP", "PP"))

# Step 4: Plot with only selected markers
plot_selected <- ggplot(data = luminex_subset_selected, aes(x = stimulation, y = log10(value))) +
  geom_violin(trim = FALSE, scale = "width") +  # Violin plots for distribution
  geom_quasirandom(size = 0.8, aes(col = HiIL6, alpha = HiIL6)) +  # Points colored by HiIL6
  facet_grid(marker ~ group, scales = "free") +  # Separate by marker and group
  scale_color_manual(values = c("green3", "red")) +  # Color scale
  scale_alpha_manual(values = c(0.5, 0.8)) +  # Transparency scale
  theme_bw() +  # Clean theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x labels
  ggtitle("Erik & Nerissa: Individuals coloured red if IL6 > 3 in the medium condition") +
  ylab("log10(concentration [pg/ml])")

# Step 5: Show plot
plot_selected


#all markers of the same ID from Erik & Nerissa
# Step 1: Set factor levels for stimulation
luminex_subset$stimulation <- factor(luminex_subset$stimulation,
                                     levels = c("M", "LPS", "SP", "PP"))

# Step 2: Plot with all markers
plot_all <- ggplot(data = luminex_subset, aes(x = stimulation, y = level)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_quasirandom(size = 0.8, aes(col = HiIL6, alpha = HiIL6)) +
  facet_grid(marker ~ group, scales = "free") +
  scale_color_manual(values = c("green3", "red")) +
  scale_alpha_manual(values = c(0.5, 0.8)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Erik & Nerissa: Individuals coloured red if IL6 > 3 in the medium condition") +
  ylab("log10(concentration [pg/ml])")

# Step 3: Show plot
plot_all




