rm(list = ls())

# install.packages("pacman")
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
               tableone, kmed, survival, survminer, skimr, lattice, reshape2, grDevices, ggrepel)

direct_plasma_patients <- read.csv("Original_data/first_syndecan_merged_markers_raw.csv")
direct_plasma_EBpatients <- direct_plasma_patients[grepl("EB", direct_plasma_patients$New_ID),]
direct_plasma_EBpatients <- direct_plasma_EBpatients %>%
  mutate(ID = str_extract(New_ID, "\\d+$"))



direct_plasma_ctr <- read.csv("Original_data/first_syndecan_list_of_healthy_raw.csv")
direct_plasma_ctr <- direct_plasma_ctr %>%
  mutate(ID = str_extract(New_ID, "\\d{4}"))


#CAP####
#### how many have whole blood stimulatioon data AND direct plasma luminex
#need to run the code what_is_wrong_long_data
table(direct_plasma_EBpatients$ID %in% luminex_long$ID )
table( unique(luminex_long$ID) %in% direct_plasma_EBpatients$ID   )


luminex_long_plamsma <- merge( direct_plasma_EBpatients, unique(luminex_long[,c("ID","HiIL6")]) , by="ID")

table(luminex_long_plamsma$HiIL6)

boxplot( log(luminex_long_plamsma$CRP) ~ luminex_long_plamsma$HiIL6 )

wilcox.test(luminex_long_plamsma$CRP ~ luminex_long_plamsma$HiIL6)

boxplot( log(luminex_long_plamsma$IL_6 ) ~ luminex_long_plamsma$HiIL6 )
wilcox.test(luminex_long_plamsma$IL_6 ~ luminex_long_plamsma$HiIL6)


boxplot( log(luminex_long_plamsma$Ferritin ) ~ luminex_long_plamsma$HiIL6 )
wilcox.test(luminex_long_plamsma$IL_6 ~ luminex_long_plamsma$HiIL6)

library(dplyr)
library(tidyr)

# Pivot biomarker columns into long format
luminex_long_biomarkers <- luminex_long_plamsma %>%
  pivot_longer(
    cols = -c(ID, New_ID, Group, HiIL6),
    names_to = "Biomarker",
    values_to = "Value"
  )

luminex_long_biomarkers$level <- log10(luminex_long_biomarkers$Value)


ggplot(data=luminex_long_biomarkers, aes(x=HiIL6, y= level))+
  geom_violin()+
  facet_wrap(~Biomarker, nrow = 4)


library(ggplot2)
library(ggpubr)

ggplot(data = luminex_long_biomarkers, aes(x = HiIL6, y = level)) +
  geom_violin(trim = FALSE) +
  facet_wrap(~Biomarker, nrow = 4, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  theme_bw()+
  theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())

#only show significant markers
library(ggpubr)

# 
p_values <- compare_means(
  formula = level ~ HiIL6,
  data = luminex_long_biomarkers,
  group.by = "Biomarker",
  method = "wilcox.test"
)

# 
p_values$p_adj <- p.adjust(p_values$p, method = "BH")

# 
sig_biomarkers <- p_values %>%
  filter(p_adj < 0.1) %>%
  arrange(p_adj)

# 
print(sig_biomarkers)



###if not BH the p values just use the raw p values##then i draw the figure
library(ggpubr)

p_values <- compare_means(
  formula = level ~ HiIL6,
  data = luminex_long_biomarkers,
  group.by = "Biomarker",
  method = "wilcox.test"
)

sig_biomarkers <- p_values %>%
  filter(p < 0.1) %>%
  arrange(p)

print(sig_biomarkers)



luminex_filtered <- luminex_long_biomarkers %>%
  filter(Biomarker %in% sig_biomarkers$Biomarker)

ggplot(data = luminex_filtered, aes(x = HiIL6, y = level)) +
  geom_violin(trim = FALSE, fill = "gray95", color = "black") +
  geom_jitter(aes(color = HiIL6), width = 0.15, size = 1.2, alpha = 0.7) +
  facet_wrap(~Biomarker, nrow = 2, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.y.npc = 0.95) +  # 👈
  scale_color_manual(values = c("LowIL6_M" = "forestgreen", "HighIL6_M" = "firebrick")) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.background = element_rect(fill = "gray90", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


#HV####
#### how many have whole blood stimulatioon data AND direct plasma luminex

table(direct_plasma_ctr$ID %in% luminex_long$ID )
table( unique(luminex_long$ID) %in% direct_plasma_ctr$ID   )


luminex_long_plamsma_ctr <- merge( direct_plasma_ctr, unique(luminex_long[,c("ID","HiIL6")]) , by="ID")

table(luminex_long_plamsma_ctr$HiIL6)

boxplot( log(luminex_long_plamsma_ctr$CRP) ~ luminex_long_plamsma_ctr$HiIL6 )

wilcox.test(luminex_long_plamsma_ctr$CRP ~ luminex_long_plamsma_ctr$HiIL6)

boxplot( log(luminex_long_plamsma_ctr$IL_6 ) ~ luminex_long_plamsma_ctr$HiIL6 )
wilcox.test(luminex_long_plamsma_ctr$IL_6 ~ luminex_long_plamsma_ctr$HiIL6)


boxplot( log(luminex_long_plamsma_ctr$Ferritin ) ~ luminex_long_plamsma_ctr$HiIL6 )
wilcox.test(luminex_long_plamsma_ctr$IL_6 ~ luminex_long_plamsma_ctr$HiIL6)

library(dplyr)
library(tidyr)

# Pivot biomarker columns into long format
luminex_long_biomarkers_ctr <- luminex_long_plamsma_ctr %>%
  pivot_longer(
    cols = -c(ID, New_ID, HiIL6),
    names_to = "Biomarker",
    values_to = "Value"
  )

luminex_long_biomarkers_ctr$level <- log10(luminex_long_biomarkers_ctr$Value)


ggplot(data=luminex_long_biomarkers_ctr, aes(x=HiIL6, y= level))+
  geom_violin()+
  facet_wrap(~Biomarker, nrow = 4)


library(ggplot2)
library(ggpubr)

ggplot(data = luminex_long_biomarkers_ctr, aes(x = HiIL6, y = level)) +
  geom_violin(trim = FALSE) +
  facet_wrap(~Biomarker, nrow = 4, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  theme_bw()+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

#only show significant markers
library(ggpubr)

library(ggpubr)

# 
p_values <- compare_means(
  formula = level ~ HiIL6,
  data = luminex_long_biomarkers_ctr,
  group.by = "Biomarker",
  method = "wilcox.test"
)

# 
p_values$p_adj <- p.adjust(p_values$p, method = "BH")

# 
sig_biomarkers <- p_values %>%
  filter(p_adj < 0.1) %>%
  arrange(p_adj)

# 
print(sig_biomarkers)


luminex_filtered_ctr <- luminex_long_biomarkers_ctr %>%
  filter(Biomarker %in% sig_biomarkers$Biomarker)

ggplot(data = luminex_filtered_ctr, aes(x = HiIL6, y = level)) +
  geom_violin(trim = FALSE, fill = "gray95", color = "black") +
  geom_jitter(aes(color = HiIL6), width = 0.15, size = 1.2, alpha = 0.7) +
  facet_wrap(~Biomarker, nrow = 2, scales = "free_y") +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.y.npc = 0.95) +  # 👈
  scale_color_manual(values = c("LowIL6_M" = "forestgreen", "HighIL6_M" = "firebrick")) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.background = element_rect(fill = "gray90", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
