rm(list = ls())

library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(grid)

# ==== 1. Load data ====
studydata_patients <- read.csv("Original_data/studydata_patients.csv")
Luminex_CAP_COVID_HV <- read.csv("Original_data/Luminex_CAP_COVID_HV.csv")
merged_data <- merge(Luminex_CAP_COVID_HV, studydata_patients, by.x = "ID", by.y = "EB_id")

# ==== 2. Long format & log10 ====
mdata <- merged_data[, 1:12]
mdata$Day <- NULL
ldata <- reshape2::melt(mdata, id = c("ID", "stimulation"))
ldata$level <- log10(ldata$value)

##export form
destination_folder <- "Original_data/" 
export_file_name <- "ldata_CAP.csv" 
write.csv(ldata, file.path(destination_folder, export_file_name), row.names = FALSE)
cat("File exported to:", file.path(destination_folder, export_file_name), "\n")




# ==== 3. Select LPS + markers ====
all_markers <- unique(ldata$variable)
pca_markers <- c("IL_1RA", "IL_1beta", "TNF", "IL_6")

ldata_sel <- ldata %>%
  filter(stimulation == "LPS", variable %in% all_markers) %>%
  filter(!is.na(level))

# ==== 4. Prepare PCA input ====
ldata_pca <- ldata_sel %>% filter(variable %in% pca_markers)
ldata_wide <- reshape2::dcast(ldata_pca, ID ~ variable, value.var = "level")
rownames(ldata_wide) <- ldata_wide$ID
ldata_wide$ID <- NULL

# ==== 5. PCA & clustering ====
pca_res <- prcomp(ldata_wide, scale. = TRUE)
kmeansres <- kmeans(pca_res$x[, 1:2], centers = 2)

# ==== 6. Format patient PCA results ====
pca_df <- data.frame(pca_res$x[, 1:2])
pca_df$ID <- as.integer(rownames(pca_df))
pca_df$cluster <- as.factor(kmeansres$cluster)
pca_df$group <- "Patient"

# ==== 7. Project HV onto PCA ====
studydata_HV <- read.csv("Original_data/studydata_HV.csv")
merged_HV <- merge(Luminex_CAP_COVID_HV, studydata_HV, by.x = "ID", by.y = "EB_id")
mdata_HV <- merged_HV[, 1:12]
mdata_HV$Day <- NULL
ldata_HV <- reshape2::melt(mdata_HV, id = c("ID", "stimulation"))
ldata_HV$level <- log10(ldata_HV$value)

##export form
destination_folder <- "Original_data/" 
export_file_name <- "ldata_HV.csv" 
write.csv(ldata_HV, file.path(destination_folder, export_file_name), row.names = FALSE)
cat("File exported to:", file.path(destination_folder, export_file_name), "\n")



ldata_HV_sel <- ldata_HV %>%
  filter(stimulation == "LPS", variable %in% pca_markers) %>%
  filter(!is.na(level))

ldata_HV_wide <- reshape2::dcast(ldata_HV_sel, ID ~ variable, value.var = "level")
rownames(ldata_HV_wide) <- ldata_HV_wide$ID
ldata_HV_wide$ID <- NULL
pca_coords_HV <- predict(pca_res, newdata = ldata_HV_wide)

pca_HV_df <- data.frame(pca_coords_HV[, 1:2])
pca_HV_df$ID <- as.integer(rownames(pca_HV_df))
pca_HV_df$group <- "Healthy group"
pca_HV_df$cluster <- "HV"

pca_all <- bind_rows(pca_df, pca_HV_df)

# ==== 8. Loadings (arrow) ====
loadings_df <- as.data.frame(pca_res$rotation[, 1:2])
loadings_df$variable <- rownames(loadings_df)
loadings_df$stim <- "LPS"
loadings_df$marker <- gsub("_", "-", loadings_df$variable)
loadings_df$label <- loadings_df$marker

# ==== 9. PCA plot ====
arrow_scale <- 6
label_offset <- 1.1
label_size <- 5
marker_colors <- c("IL-1RA" = "#4DAF4A", "IL-1beta" = "#FF6347", "TNF" = "#FF7F00", "IL-6" = "#F781BF")
stim_cols_arrow <- c("LPS" = "#FF7F00")

ggplot(pca_all, aes(PC1, PC2, color = cluster, shape = group)) +
  geom_point(size = 2, alpha = 0.85) +
  stat_ellipse(aes(group = cluster), linewidth = 1.2) +
  geom_segment(data = loadings_df,
               aes(x = 0, y = 0,
                   xend = PC1 * arrow_scale,
                   yend = PC2 * arrow_scale,
                   color = stim),
               arrow = arrow(length = unit(0.35, "cm")),
               linewidth = 1.2,
               inherit.aes = FALSE) +
  geom_label_repel(data = loadings_df,
                   aes(x = PC1 * arrow_scale * label_offset,
                       y = PC2 * arrow_scale * label_offset,
                       label = label, fill = marker),
                   color = "black",
                   size = label_size,
                   label.padding = 0.3,
                   label.size = 0.4,
                   segment.color = "grey30",
                   max.overlaps = Inf,
                   box.padding = 0.5,
                   inherit.aes = FALSE) +
  scale_color_manual(
    values = c("1" = "#E41A1C", "2" = "#4DAF4A", "HV" = "#377EB8", stim_cols_arrow),
    labels = c("Cluster 1", "Cluster 2", "Healthy group", "LPS stimulus")
  ) +
  scale_fill_manual(values = marker_colors, name = "Marker") +
  scale_shape_manual(values = c("Patient" = 17, "Healthy group" = 16)) +
  labs(
    title = "PCA projection based on LPS: IL-1RA, IL-1β, TNF, IL-6",
    x = paste0("PC1 (", round(summary(pca_res)$importance[2, 1] * 100, 2), "%)"),
    y = paste0("PC2 (", round(summary(pca_res)$importance[2, 2] * 100, 2), "%)"),
    color = "Cluster / Stimulus",
    shape = "Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# ==== 10. Boxplot of all marker expression (add HV) ====
cluster_df <- pca_df[, c("ID", "cluster")]
cluster_df$ID <- as.integer(cluster_df$ID)

ldata_all_box_patients <- ldata_sel %>%
  filter(variable %in% all_markers) %>%
  inner_join(cluster_df, by = "ID")

ldata_all_box_HV <- ldata_HV %>%
  filter(stimulation == "LPS", variable %in% all_markers) %>%
  filter(!is.na(level)) %>%
  mutate(cluster = "HV")

ldata_all_box_combined <- bind_rows(ldata_all_box_patients, ldata_all_box_HV)

ggplot(ldata_all_box_combined, aes(x = cluster, y = level, fill = cluster)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75)) +
  geom_jitter(shape = 21, color = "black", alpha = 0.3,
              position = position_dodge(width = 0.75), stroke = 0.3) +
  facet_wrap(~ variable, ncol = 4, scales = "free_y") +
  scale_fill_manual(values = c("HV" = "#377EB8", "1" = "#E41A1C", "2" = "#4DAF4A")) +
  theme_bw(base_size = 14) +
  labs(
    title = "All Marker Expression by PCA Cluster (LPS)",
    x = "Cluster", y = "Log10 Marker Level"
  ) +
  theme(
    strip.text = element_text(size = 12),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    axis.title = element_text(size = 13),
    legend.position = "none"
  )

# ==== 11. Boxplot of PC1 and PC2 ====
pca_all$Group3 <- factor(pca_all$cluster, levels = c("HV", "1", "2"),
                         labels = c("Healthy", "Cluster 1", "Cluster 2"))

pca_all_long <- pca_all %>%
  pivot_longer(cols = c(PC1, PC2), names_to = "Component", values_to = "Value")

ggplot(pca_all_long, aes(x = Group3, y = Value, fill = Group3)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4) +
  facet_wrap(~ Component, scales = "free_y") +
  scale_fill_manual(values = c("Healthy" = "#377EB8", "Cluster 1" = "#E41A1C", "Cluster 2" = "#4DAF4A")) +
  labs(
    title = "PC1 and PC2 Distribution across Groups",
    x = "Group", y = "PCA Value"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(size = 13),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    legend.position = "none"
  )


#marker the p value
library(ggpubr)  # 用于添加 p-value

# ==== 10. Boxplot of all marker expression with p-value ====
cluster_df <- pca_df[, c("ID", "cluster")]
cluster_df$ID <- as.integer(cluster_df$ID)

ldata_all_box_patients <- ldata_sel %>%
  filter(variable %in% all_markers) %>%
  inner_join(cluster_df, by = "ID")

ldata_all_box_HV <- ldata_HV %>%
  filter(stimulation == "LPS", variable %in% all_markers) %>%
  filter(!is.na(level)) %>%
  mutate(cluster = "HV")

ldata_all_box_combined <- bind_rows(ldata_all_box_patients, ldata_all_box_HV)

# 指定比较组
my_comparisons <- list(c("HV", "1"), c("HV", "2"), c("1", "2"))

# 画图并加统计显著性
ggplot(ldata_all_box_combined, aes(x = cluster, y = level, fill = cluster)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75)) +
  geom_jitter(shape = 21, color = "black", alpha = 0.3,
              position = position_dodge(width = 0.75), stroke = 0.3) +
  facet_wrap(~ variable, ncol = 4, scales = "free_y") +
  scale_fill_manual(values = c("HV" = "#377EB8", "1" = "#E41A1C", "2" = "#4DAF4A")) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", 
                     label = "p.signif", step.increase = 0.08, size = 3) +
  theme_bw(base_size = 14) +
  labs(
    title = "All Marker Expression by PCA Cluster (LPS)",
    y = "Log10 Marker Level"  
  ) +
  theme(
    strip.text = element_text(size = 12),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    axis.title.x = element_blank(),  
    axis.title.y = element_text(size = 13),
    legend.position = "none"
  )

# 
ggsave("Figure/Boxplot_Marker_by_Cluster_LPS.svg",
       width = 14, height = 8, units = "in", dpi = 300)


library(ggpubr)

my_comparisons_pca <- list(c("Healthy", "Cluster 1"),
                           c("Healthy", "Cluster 2"),
                           c("Cluster 1", "Cluster 2"))

ggplot(pca_all_long, aes(x = Group3, y = Value, fill = Group3)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4) +
  facet_wrap(~ Component, scales = "free_y") +
  scale_fill_manual(values = c("Healthy" = "#377EB8", "Cluster 1" = "#E41A1C", "Cluster 2" = "#4DAF4A")) +
  stat_compare_means(
    comparisons = my_comparisons_pca,
    method = "wilcox.test",
    label = "p.signif",
    step.increase = 0.08,
    size = 3,
    tip.length = 0.02,
    symnum.args = list(
      cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
      symbols = c("***", "***", "**", "*", "ns")  # 限制最多三个星号
    )
  ) +
  labs(
    title = "PC1 and PC2 Distribution across Groups",
    x = "Group", y = "PCA Value"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(size = 13),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    legend.position = "none",
    plot.margin = margin(10, 20, 30, 10)
  ) +
  coord_cartesian(clip = "off")



plot(ldata_wide$TNF, ldata_wide$IL_1RA)
points(ldata_HV_wide$TNF, ldata_HV_wide$IL_1RA, pch=16, col="green")


t.test(ldata_wide$TNF , ldata_HV_wide$TNF )


ldata_HV_sel <- ldata_HV[ldata_HV$stimulation=="LPS" & ldata_HV$variable %in% c("TNF","IL_1RA"),]

ldata_HV_wide <- reshape2::dcast(ldata_HV_sel, ID ~ variable, value.var = "level")


points(ldata_HV_wide$TNF, ldata_HV_wide$IL_1RA, pch=16, col="green")
plot(ldata_HV_wide$TNF, ldata_HV_wide$IL_1RA, pch=16, col="green")
cor(ldata_HV_wide$TNF, ldata_HV_wide$IL_1RA, pch=16, col="green")
cor(ldata_HV_wide$TNF, ldata_HV_wide$IL_1RA)


ldata_wide$hi_IL1RA <- ifelse(ldata_wide$IL_1RA > 4.318, "hiIL1RA", "lowILRA")
ldata_wide$hi_TNF <- ifelse(ldata_wide$TNF < 2.9024, "lowTNF", "hiTNF")

table(ldata_wide$hi_IL1RA , ldata_wide$hi_TNF)

ldata_wide$group2 <- paste0(ldata_wide$hi_IL1RA, ldata_wide$hi_TNF )

plot(ldata_wide$TNF, ldata_wide$IL_1RA, 
     col=as.factor(ldata_wide$group2), pch=16)






