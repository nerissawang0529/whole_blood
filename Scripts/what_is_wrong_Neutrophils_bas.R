# rm(list = ls())
# need to run the codes about what is wrong first

# install.packages("pacman")
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
               tableone, kmed, survival, survminer, skimr, lattice, reshape2, grDevices, ggrepel)

Neutrophils <- readxl::read_xlsx("Original_data/Bas_combined_PMN_stimulation.xlsx")
library(readxl)
library(dplyr)

# Assign full data as main_df
main_df <- Neutrophils

# Load redo measurements
redo1 <- read_excel("Original_data/Bas_combined_PMN_stimulation.xlsx", sheet = "Redo samples")
redo2 <- read_excel("Original_data/Bas_combined_PMN_stimulation.xlsx", sheet = "Redo samples rest")

# Convert all values to character to avoid type mismatches
redo1[] <- lapply(redo1, as.character)
redo2[] <- lapply(redo2, as.character)

# Combine redo datasets
redo_df <- bind_rows(redo1, redo2)

# Define unique key columns
key_cols <- c("ID", "Inhoud")

# Select marker names
biomarkers <- c("IL8", "Lipocalin2_NGAL", "Proteinase3", "Lactoferrin", "MPO")

# Replace OOR values in main_df with redo values if available and valid
for (marker in biomarkers) {
  oor_rows <- which(main_df[[marker]] == "OOR >")
  
  for (i in oor_rows) {
    id <- main_df$ID[i]
    inhoud <- main_df$Inhoud[i]
    
    replacement <- redo_df %>%
      filter(ID == id, Inhoud == inhoud) %>%
      pull(marker)
    
    if (length(replacement) > 0 && !grepl("OOR", replacement[1])) {
      main_df[[marker]][i] <- replacement[1]
    }
  }
}

# Exclude D28 samples
main_df <- main_df[!grepl("_D28$", main_df$ID), ]

# Merge with luminex_long
main_df$ID <- as.character(main_df$ID)
luminex_long$ID <- as.character(luminex_long$ID)

# Select metadata info (group and IL-6 subgroup)
meta_info <- luminex_long %>% select(ID, group, HiIL6) %>% distinct()

# Merge into main_df and filter missing group values
main_df <- left_join(main_df, meta_info, by = "ID")
main_df <- main_df %>% filter(!is.na(group))

# Ensure factor order for IL-6 grouping
main_df$HiIL6 <- factor(main_df$HiIL6, levels = c("LowIL6_M", "HighIL6_M"))

##need to change the group##############################################################################
main_df <- main_df %>% filter(group == "HV")


# Filter only M_sup condition (unstimulated) for violin plots
m_data <- main_df %>%
  filter(grepl("M_sup", Inhoud)) %>%
  select(ID, HiIL6, IL8, Lipocalin2_NGAL, Proteinase3, Lactoferrin, MPO)

# Long format + log10 transform
m_long <- m_data %>%
  pivot_longer(cols = -c(ID, HiIL6), names_to = "marker", values_to = "value") %>%
  dplyr::mutate(value = as.numeric(value)) %>%
  dplyr::filter(!is.na(value), value > 0) %>%
  dplyr::mutate(log_value = log10(value))

m_long <- m_long %>%
  filter(marker %in% c("Lipocalin2_NGAL", "Proteinase3", "MPO"))


# Violin plots of unstimulated values by IL-6 group
ggplot(m_long, aes(x = HiIL6, y = log_value, fill = HiIL6)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  facet_wrap(~marker, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  scale_fill_manual(values = c("HighIL6_M" = "#E74C3C",  # Red
                               "LowIL6_M"  = "#3498DB")) +  # Blue
  theme_minimal() +
  labs(y = "log10(Value)", x = "Unstimulated IL-6 group",
       title = "Neutrophil degranulation markers under unstimulated conditions in IL-6 groups defined by whole blood IL-6 levels in medium condition (non-infection)")








# === Δ release calculation section ===

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

# 1. Select neutrophil data under M / LPS / KP / SP stimulation
delta_data <- main_df %>%
  filter(grepl("M_sup|LPS_sup|KP_sup|SP_sup", Inhoud)) %>%
  select(ID, Inhoud, HiIL6, IL8, Lipocalin2_NGAL, Proteinase3, Lactoferrin, MPO)

# 2. Simplify Inhoud labels
delta_data <- delta_data %>%
  mutate(Inhoud_clean = case_when(
    grepl("M_sup", Inhoud)   ~ "M",
    grepl("LPS_sup", Inhoud) ~ "LPS",
    grepl("KP_sup", Inhoud)  ~ "KP",
    grepl("SP_sup", Inhoud)  ~ "SP",
    TRUE ~ "Other"
  ))

# 3. Reshape → log10 → pivot wide
delta_long <- delta_data %>%
  pivot_longer(cols = -c(ID, HiIL6, Inhoud, Inhoud_clean),
               names_to = "marker", values_to = "value") %>%
  mutate(value = as.numeric(value)) %>%
  filter(!is.na(value), value > 0) %>%
  mutate(log_value = log10(value)) %>%
  select(ID, HiIL6, marker, Inhoud_clean, log_value) %>%
  distinct(ID, HiIL6, marker, Inhoud_clean, .keep_all = TRUE) %>%
  pivot_wider(names_from = Inhoud_clean, values_from = log_value)

# 4. Compute Δ log10(stim - M)
delta_long <- delta_long %>%
  mutate(
    delta_LPS = as.numeric(LPS) - as.numeric(M),
    delta_KP  = as.numeric(KP)  - as.numeric(M),
    delta_SP  = as.numeric(SP)  - as.numeric(M)
  )

# 5. Prepare for plotting
plot_data <- delta_long %>%
  pivot_longer(cols = c(delta_LPS, delta_KP, delta_SP),
               names_to = "stimulus", values_to = "delta_log") %>%
  filter(!is.na(delta_log))

# 6. Violin plots of Δ stimulated - baseline (M)
ggplot(plot_data, aes(x = HiIL6, y = delta_log, fill = HiIL6)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  facet_wrap(~stimulus + marker, scales = "free_y", ncol = 5) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  scale_fill_manual(values = c("HighIL6_M" = "#E74C3C", "LowIL6_M" = "#3498DB")) +
  theme_minimal(base_size = 14) +
  labs(
    y = expression(Delta~log[10]~"(stimulated - M_sup)"),
    x = "IL-6 group (based on whole blood IL-6 in medium condition)",
    title = "Relative increase in neutrophil degranulation markers upon stimulation, normalized to baseline"
  )
