#need to run what_is_wrong_long_data, and what_is_wrong first to get the database

### high. versus low IL6

highIL6 <- datab[log10(datab$IL_6) > 3 & datab$stimulation=="M" , ]

table(highIL6$group)

table(datab$group[ datab$stimulation=="M" ])

45/221
13/175
27/93

luminex_long$HiIL6 <- ifelse( luminex_long$ID %in% highIL6$ID, 
                              "HighIL6_M", "LowIL6_M")



table(luminex_long$HiIL6)

3060/17602


luminex_long$value <- as.numeric(luminex_long$value)
luminex_long$level <- log10(luminex_long$value)
luminex_long$stimulation <- as.factor(luminex_long$stimulation) 
luminex_long$stimulation <- factor(luminex_long$stimulation, levels=c("M","LPS","KP","PP"))

luminex_long$group<- as.factor(luminex_long$group) 
luminex_long$group<- factor(luminex_long$group, levels=c("HV","CAP","COVID"))

library(ggbeeswarm)

plot <- ggplot( data=luminex_long, aes(x=stimulation, y=level) )+
  geom_violin(trim = FALSE, scale = "width")+
  geom_quasirandom(size=0.8,  aes(col=HiIL6, alpha=HiIL6))+
  facet_grid(marker ~ group, scales="free")+
  scale_color_manual(values=c("green3","red"))+
  scale_alpha_manual(values=c(0.5, 0.8))+
  theme_bw()+
  ggtitle("Individuals coloured red if IL6 > 3 in the medium condition ")+
  ylab("log10(concentration[pg/ml])")

plot
ggsave(file="/Users/huiwang/Documents/Luminex/luminex_R/Figure/groupVmarker1.svg", plot=plot, width=11, height=14)


luminex_long$HiIL6 <- as.factor(luminex_long$HiIL6) 
luminex_long$HiIL6<- factor(luminex_long$HiIL6, levels=c("LowIL6_M","HighIL6_M"))


plot2 <- ggplot( data=luminex_long, aes(x= stimulation, y=level) )+
  #geom_violin(trim = FALSE, scale = "width")+
  geom_point( aes(group=ID) )+
  facet_grid(marker ~ group*HiIL6, scales="free")+
  #scale_color_manual(values=c("green3","red"))+
  #scale_alpha_manual(values=c(0.5, 0.8))+
  theme_bw()+
  ggtitle("Individuals coloured red if IL6 > 3 in the medium condition ")+
  ylab("log10(concentration[pg/ml])")

plot2

plot2 <- ggplot(data = luminex_long, aes(x = stimulation, y = level, group = ID)) +
  geom_line(alpha = 0.5, color = "gray") +  # Connect dots per ID
  geom_point(aes(color = HiIL6), size = 2) +  # Color by IL-6 high/low
  facet_grid(marker ~ group * HiIL6, scales = "free") +
  theme_bw() +
  ggtitle("Markers levels per stimulation - Lines connect same individual (ID)") +
  ylab("log10(concentration [pg/ml])") 
  
plot2


ggsave(file="/Users/huiwang/Documents/Luminex/luminex_R/Figure/groupVmarker2.svg", plot=plot2, width=11, height=14)

#significant difference for plot2
# Create an empty list to store results
unique_groups <- unique(luminex_long$group)
unique_il6 <- unique(luminex_long$HiIL6)
unique_markers <- unique(luminex_long$marker)


wilcox_results <- list()

# Start index
idx <- 1

# Loop over group, HiIL6, and marker
for (grp in unique(luminex_long$group)) {
  for (il6 in unique(luminex_long$HiIL6)) {
    for (mk in unique(luminex_long$marker)) {
      
      subset_df <- luminex_long[
        luminex_long$group == grp &
          luminex_long$HiIL6 == il6 &
          luminex_long$marker == mk, 
      ]
      
      if (length(unique(subset_df$stimulation)) > 1) {
        test <- pairwise.wilcox.test(
          x = subset_df$level,
          g = subset_df$stimulation,
          p.adjust.method = "bonferroni"
        )
        
        # Extract pairwise results
        pw_matrix <- as.data.frame(as.table(test$p.value))
        colnames(pw_matrix) <- c("Group1", "Group2", "p_value")
        pw_matrix <- na.omit(pw_matrix)
        pw_matrix$Group <- grp
        pw_matrix$HiIL6 <- il6
        pw_matrix$Marker <- mk
        
        # Store in list
        wilcox_results[[idx]] <- pw_matrix
        idx <- idx + 1
      }
    }
  }
}

# Combine into one data frame
final_results <- do.call(rbind, wilcox_results)

# Optional: Order columns for clarity
final_results <- final_results[, c("Group", "HiIL6", "Marker", "Group1", "Group2", "p_value")]

# View first few lines
head(final_results)


#only for HV AND CAP groups of the p values ####
# Remove rows where group is "COVID"
luminex_long_1 <- luminex_long[luminex_long$group != "COVID", ]
plot21 <- ggplot(data = luminex_long_1, aes(x = stimulation, y = level, group = ID)) +
  geom_line(alpha = 0.4, color = "gray70") +  # lighter connecting lines
  geom_point(aes(color = HiIL6), size = 2) +  # colored points
  facet_grid(marker ~ group * HiIL6, scales = "free") +
  scale_color_manual(
    values = c("LowIL6_M" = "#377EB8",  # blue (colorblind-friendly)
               "HighIL6_M" = "#E41A1C") # vermillion/red
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0, face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  ) +
  ggtitle("Markers levels per stimulation - Lines connect same individual (ID)") +
  ylab("log10(concentration [pg/ml])") 

plot21
ggsave(file="/Users/huiwang/Documents/Luminex/luminex_R/Figure/group_stimulation.svg", plot=plot21, width=12, height=12)

#export 
destination_folder <- "Original_data/" 
export_file_name <- "luminex_long_1.csv" 
write.csv(luminex_long_1, file.path(destination_folder, export_file_name), row.names = FALSE)
cat("File exported to:", file.path(destination_folder, export_file_name), "\n")

#significant difference for plot2
# Create an empty list to store results
unique_groups <- unique(luminex_long_1$group)
unique_il6 <- unique(luminex_long_1$HiIL6)
unique_markers <- unique(luminex_long_1$marker)


wilcox_results <- list()

# Start index
idx <- 1

# Loop over group, HiIL6, and marker
for (grp in unique(luminex_long_1$group)) {
  for (il6 in unique(luminex_long_1$HiIL6)) {
    for (mk in unique(luminex_long_1$marker)) {
      
      subset_df <- luminex_long_1[
        luminex_long_1$group == grp &
          luminex_long_1$HiIL6 == il6 &
          luminex_long_1$marker == mk, 
      ]
      
      if (length(unique(subset_df$stimulation)) > 1) {
        test <- pairwise.wilcox.test(
          x = subset_df$level,
          g = subset_df$stimulation,
          p.adjust.method = "bonferroni"
        )
        
        # Extract pairwise results
        pw_matrix <- as.data.frame(as.table(test$p.value))
        colnames(pw_matrix) <- c("Group1", "Group2", "p_value")
        pw_matrix <- na.omit(pw_matrix)
        pw_matrix$Group <- grp
        pw_matrix$HiIL6 <- il6
        pw_matrix$Marker <- mk
        
        # Store in list
        wilcox_results[[idx]] <- pw_matrix
        idx <- idx + 1
      }
    }
  }
}

# Combine into one data frame
final_results_1 <- do.call(rbind, wilcox_results)

# Optional: Order columns for clarity
final_results_1 <- final_results_1[, c("Group", "HiIL6", "Marker", "Group1", "Group2", "p_value")]
# Add a column to flag significance
final_results_1$significant <- final_results_1$p_value < 0.05

#export 
destination_folder <- "Original_data/" 
export_file_name <- "p_patientsgroup_stimulation.csv" 
write.csv(final_results_1, file.path(destination_folder, export_file_name), row.names = FALSE)
cat("File exported to:", file.path(destination_folder, export_file_name), "\n")


#####




### correlation facets
library(dplyr)
library(tidyr)

marker_stim_wide <- luminex_long %>%
  dplyr::select(ID, marker, stimulation, level, group, HiIL6) %>%
  tidyr::pivot_wider(
    names_from = stimulation,
    values_from = level
  )

head(marker_stim_wide)



marker_stim_long <- marker_stim_wide %>%
  pivot_longer(
    cols = c(KP, LPS, PP),              # the stim columns to melt
    names_to = "stimulation",           # new categorical column
    values_to = "stimulation_values"    # new value column
  )

head(marker_stim_long)


plot3 <- ggplot(data = marker_stim_long, aes(x = M, y = stimulation_values)) +
  geom_point( size = 1, aes(col=HiIL6)) +  # Color by IL-6 high/low
  facet_grid(marker ~ group * stimulation, scales = "free") +
  theme_bw() +
  ggtitle("IL-8 levels per stimulation - Lines connect same individual (ID)") +
  ylab("log10(concentration [pg/ml])") +
  scale_color_manual(values=c("green3","red"))

plot3


ggsave(file="/Users/huiwang/Documents/Luminex/luminex_R/Figure/correlations2.svg", plot=plot3, width=15, height=20)



plot3 <- ggplot(data = marker_stim_long, aes(x = M, y = stimulation_values)) +
  geom_point( size = 1, aes(col=M>2)) +  # Color by IL-6 high/low
  facet_grid(marker ~ group * stimulation, scales = "free") +
  theme_bw() +
  ggtitle("IL-8 levels per stimulation - Lines connect same individual (ID)") +
  ylab("log10(concentration [pg/ml])") 
  #scale_color_manual(values=c("green3","red"))

plot3


ggsave(file="/Users/huiwang/Documents/Luminex/luminex_R/Figure/correlations2.svg", plot=plot3, width=15, height=20)

### delta

marker_stim_long$delta <- marker_stim_long$stimulation_values - marker_stim_long$M

hist(marker_stim_long$delta)

plot4 <- ggplot(data = marker_stim_long, aes(x = stimulation, y = delta)) +
  geom_violin(trim = FALSE, scale = "width")+
  geom_quasirandom(size=0.8,  aes(col=HiIL6))+
  facet_grid(marker ~ group , scales = "free") +
  theme_bw() +
  ggtitle("IL-8 levels per stimulation - Lines connect same individual (ID)") +
  ylab("log10(concentration [pg/ml])") 
#scale_color_manual(values=c("green3","red"))

plot4

ggsave(file="/Users/huiwang/Documents/Luminex/luminex_R/Figure/delta.svg", plot=plot4, width=15, height=20)

#compare CAP AND HV in m
# 
install.packages(c("dplyr", "ggpubr", "rstatix"))

# 
library(dplyr)
library(ggpubr)
library(rstatix)

# 1. 
df_M <- luminex_long %>%
  filter(stimulation == "M", group %in% c("CAP", "HV")) %>%
  mutate(group = recode(group,
                        "HV" = "Non-infection",
                        "CAP" = "CAP"))  

# 2. 
stat_test <- df_M %>%
  group_by(marker) %>%
  wilcox_test(level ~ group) %>%
  adjust_pvalue(method = "BH") %>%
  mutate(p.adj.signif = case_when(
    p.adj <= 0.001 ~ "***",
    p.adj <= 0.01 ~ "**",
    p.adj <= 0.05 ~ "*",
    TRUE ~ "ns"
  ))

# 
stat_test <- stat_test %>% add_xy_position(x = "group")

# 
ggviolin(df_M, x = "group", y = "level", fill = "group",
         palette = "jco", add = "boxplot", facet.by = "marker",
         ylab = "Cytokine value", xlab = "Group") +
  stat_pvalue_manual(stat_test, label = "p.adj.signif",
                     tip.length = 0.01, hide.ns = TRUE) +
  theme(strip.text = element_text(size = 12, face = "bold"))


