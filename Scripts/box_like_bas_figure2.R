rm(list = ls())

library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(grid)
## install.packages("pacman")
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
               tableone, kmed, survival, survminer, skimr, lattice, reshape2, grDevices, ggrepel)

# ==== 1. 导入数据 ====
studydata_patients <- read.csv("Original_data/studydata_patients.csv")
Luminex_CAP_COVID_HV <- read.csv("Original_data/Luminex_CAP_COVID_HV.csv")
merged_data <- merge(Luminex_CAP_COVID_HV, studydata_patients, 
                     by.x = "ID", by.y = "EB_id")

# ==== 2. 长格式 & 取 log ====
mdata <- merged_data[, 1:12]
mdata$Day <- NULL
ldata <- reshape2::melt(mdata, id = c("ID", "stimulation"))
ldata$level <- log10(ldata$value)

library(tidyverse)

# Step 1: Select TNF under LPS stimulation and rank
tnf_lps <- ldata %>%
  filter(stimulation == "LPS", variable == "TNF") %>%
  arrange(value)

# Step 2: Define low/high TNF groups (25% lowest and 25% highest)
n <- nrow(tnf_lps)
tnf_lps <- tnf_lps %>%
  mutate(TNF_group = case_when(
    row_number() <= floor(0.25 * n) ~ "Low TNF",
    row_number() > ceiling(0.75 * n) ~ "High TNF",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(TNF_group)) %>%
  select(ID, TNF_group)

# Step 3: Merge back group info
ldata_with_group <- ldata %>%
  filter(stimulation %in% c("LPS", "KP", "PP")) %>%
  inner_join(tnf_lps, by = "ID")

# Step 4: Select cytokines
cytokines_of_interest <- c("IL_1beta", "IL_6", "IL_10", "IL_1RA")
plot_data <- ldata_with_group %>%
  filter(variable %in% cytokines_of_interest)

# Step 5: Create combined group (e.g. "LPS_High")
plot_data <- plot_data %>%
  mutate(stim_group = paste(stimulation, TNF_group, sep = "_")) %>%
  mutate(stim_group = factor(stim_group,
                             levels = c("LPS_Low TNF", "LPS_High TNF",
                                        "KP_Low TNF", "KP_High TNF",
                                        "PP_Low TNF", "PP_High TNF")),
         variable = factor(variable, levels = cytokines_of_interest))

# Step 6: Plot
library(ggpubr)

library(ggpubr)
library(ggplot2)

# Update stimulation labels: KP -> Kpneu, PP -> Pneu

library(ggpubr)
library(ggplot2)
library(dplyr)

# Ensure TNF group order
plot_data$TNF_group <- factor(plot_data$TNF_group, levels = c("Low TNF", "High TNF"))

# Comparison list
my_comparisons <- list(c("Low TNF", "High TNF"))

# Plot
p_list <- lapply(cytokines_of_interest, function(marker) {
  ggplot(filter(plot_data, variable == marker),
         aes(x = TNF_group, y = value, color = TNF_group)) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
    stat_summary(fun = median, geom = "crossbar", width = 0.5, fatten = 1.5, color = "black") +
    stat_compare_means(
      method = "wilcox.test",
      comparisons = my_comparisons,
      aes(label = ..p.signif..),
      group.by = "stimulation",
      hide.ns = TRUE
    ) +
    scale_y_log10() +
    scale_color_manual(values = c("Low TNF" = "green", "High TNF" = "red")) +
    facet_grid(~stimulation) +
    labs(title = marker, x = NULL, y = NULL, color = "TNF-α group") +
    theme_bw(base_size = 14) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.text = element_text(size = 12, face = "bold"),
      legend.position = "bottom"
    )
})

# Combine
final_plot <- cowplot::plot_grid(plotlist = p_list, ncol = 2)
print(final_plot)
