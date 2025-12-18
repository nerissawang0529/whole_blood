rm(list = ls ())

## install.packages("pacman")
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
               tableone, kmed, survival, survminer, skimr, lattice, reshape2, grDevices, ggrepel)


#ready data
studydata_HV <- read.csv("Original_data/studydata_HV.csv")

#healthy_like_patients
healthy_like_patients <- studydata_HV[
  studydata_HV$hypertension == "No" &
    studydata_HV$COPD == "No" &
    studydata_HV$diabetes == 0 &
    studydata_HV$ccd == "No" &
    studydata_HV$ckd == "No" &
    studydata_HV$mneoplasm == "No" &
    studydata_HV$immune_sup == "No" &
    studydata_HV$cnd == "No" &
    studydata_HV$cpd == "No",   
]

##export form
destination_folder <- "Original_data/" 
export_file_name <- "studydata_healthy_like_patients.csv" 
write.csv(healthy_like_patients, file.path(destination_folder, export_file_name), row.names = FALSE)
cat("File exported to:", file.path(destination_folder, export_file_name), "\n")


#run what_is_wrong get the figure IL-6 levels in non-infectin group under medium over time
# 
healthy_ids <- read.csv("Original_data/studydata_healthy_like_patients.csv")$EB_id

# 
databM_HV <- databM_HV %>%
  mutate(is_healthy = ifelse(ID %in% healthy_ids, "Healthy-like in non-infection control", "Non-Healthy in non-infection control"))

# 
ggplot(data = databM_HV, aes(x = date_created2, y = log10(IL_6), color = is_healthy)) +
  geom_point(alpha = 0.8, size = 2.2) +
  scale_color_manual(values = c("Healthy-like in non-infection control" = "black", "Non-Healthy in non-infection control" = "grey")) +
  scale_x_date(date_breaks = "3 month", date_labels = "%Y-%m") +
  labs(
    x = "Sample Collection Date (by 3 Months)",
    y = "log10(IL-6)",
    color = "Group",
    title = "IL-6 Over Time in Healthy-like vs. Non-Healthy Individuals Among Non-Infection Controls (Medium)"
  ) +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




# different markers
healthy_ids <- read.csv("Original_data/studydata_healthy_like_patients.csv")$EB_id
databM_HV <- databM %>%
  filter(group == "HV") %>%
  mutate(HealthGroup = ifelse(ID %in% healthy_ids, 
                              "Healthy-like in non-infection control", 
                              "Non-Healthy in non-infection control"))

# 
markers <- c("IL_8_1", "IL_10", "TNF", "CCL2", "CCL3", "CCL4", "IL_1RA", "IL_1beta")

# 
datab_long1 <- databM_HV %>%
  pivot_longer(cols = all_of(markers), names_to = "marker", values_to = "value") %>%
  mutate(log_value = log10(value))

# 
datab_long1$HealthGroup <- recode(datab_long1$HealthGroup,
                                  "Healthy-like in non-infection control" = "Healthy-like",
                                  "Non-Healthy in non-infection control" = "Non-Healthy")

pp <- ggplot(datab_long1, aes(x = HealthGroup, y = log_value, fill = HealthGroup)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
  geom_jitter(width = 0.15, size = 1.3, alpha = 0.4, color = "black") +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 3.5) +
  facet_wrap(~ marker, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("Healthy-like" = "black", 
                               "Non-Healthy" = "grey")) +
  theme_classic(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    x = "",
    y = "log10(Marker Value)",
    title = "Baseline Cytokine Levels in Healthy-like vs Non-Healthy Individuals (Non-Infection Controls)"
  )

pp
ggsave(file="/Users/huiwang/Documents/Luminex/luminex_R/Figure/healthy_non_healthy.svg", plot=pp, width=12, height=8)
