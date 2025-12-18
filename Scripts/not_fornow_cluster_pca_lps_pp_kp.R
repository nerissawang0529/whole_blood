rm(list = ls())

library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(grid)

# ==== 1. 导入数据 ====
studydata_patients <- read.csv("Original_data/studydata_patients.csv")
Luminex_CAP_COVID_HV <- read.csv("Original_data/Luminex_CAP_COVID_HV.csv")
merged_data <- merge(Luminex_CAP_COVID_HV, studydata_patients, 
                     by.x = "ID", by.y = "EB_id")

# ==== 2. 长格式 & 取 log ====
mdata <- merged_data[, 1:12]
mdata$Day <- NULL
ldata <- melt(mdata, id = c("ID", "stimulation"))
ldata$level <- log10(ldata$value)

# ==== 3. 筛选 LPS + 四个marker ====
selected_markers <- c("IL_1RA", "IL_1beta", "TNF", "IL_6")

ldata_sel <- ldata %>%
  filter(stimulation == "LPS", variable %in% selected_markers) %>%
  filter(!is.na(level))

# 宽格式
ldata_wide <- dcast(ldata_sel, ID ~ variable, value.var = "level")
rownames(ldata_wide) <- ldata_wide$ID
ldata_wide$ID <- NULL

# ==== 4. PCA & kmeans聚类 ====
pca_res <- prcomp(ldata_wide, scale. = TRUE)
kmeansres <- kmeans(pca_res$x[, 1:2], centers = 2)

# ==== 5. 整合cluster、PCA坐标、group信息 ====
pca_df <- data.frame(pca_res$x[,1:2])
pca_df$ID <- rownames(pca_df)
pca_df$cluster <- as.factor(kmeansres$cluster)
pca_df$group <- "Patient"

# 加入HV（只投影，不参与PCA计算）
studydata_HV <- read.csv("Original_data/studydata_HV.csv")
merged_HV <- merge(Luminex_CAP_COVID_HV, studydata_HV, by.x = "ID", by.y = "EB_id")
mdata_HV <- merged_HV[, 1:12]
mdata_HV$Day <- NULL
ldata_HV <- melt(mdata_HV, id = c("ID", "stimulation"))
ldata_HV$level <- log10(ldata_HV$value)

ldata_HV_sel <- ldata_HV %>%
  filter(stimulation == "LPS", variable %in% selected_markers) %>%
  filter(!is.na(level))

ldata_HV_wide <- dcast(ldata_HV_sel, ID ~ variable, value.var = "level")
rownames(ldata_HV_wide) <- ldata_HV_wide$ID
ldata_HV_wide$ID <- NULL
pca_coords_HV <- predict(pca_res, newdata = ldata_HV_wide)

# 合并HV和病人
pca_HV_df <- data.frame(pca_coords_HV[,1:2])
pca_HV_df$ID <- rownames(pca_HV_df)
pca_HV_df$group <- "Healthy group"
pca_HV_df$cluster <- "HV"

pca_all <- bind_rows(pca_df, pca_HV_df)

# ==== 6. 制作载荷箭头（loadings） ====
loadings_df <- as.data.frame(pca_res$rotation[, 1:2])
loadings_df$variable <- rownames(loadings_df)
loadings_df$stim <- "LPS"  # 因为只用了LPS
loadings_df$marker <- gsub("_", "-", loadings_df$variable)
loadings_df$label <- loadings_df$marker

# ==== 7. 自定义颜色 ====
marker_colors <- c(
  "IL-1RA" = "#4DAF4A",
  "IL-1beta" = "#FF6347",
  "TNF" = "#FF7F00",
  "IL-6" = "#F781BF"
)

stim_cols_arrow <- c("LPS" = "#FF7F00")

# ==== 8. 绘图 ====
arrow_scale <- 6
label_offset <- 1.1
label_size <- 5

ggplot(pca_all, aes(PC1, PC2, color = cluster, shape = group)) +
  geom_point(size = 2, alpha = 0.85) +
  stat_ellipse(aes(group = cluster), linewidth = 1.2) +
  
  # 载荷箭头
  geom_segment(data = loadings_df,
               aes(x = 0, y = 0,
                   xend = PC1 * arrow_scale,
                   yend = PC2 * arrow_scale,
                   color = stim),
               arrow = arrow(length = unit(0.35, "cm")),
               linewidth = 1.2,
               inherit.aes = FALSE) +
  
  # 标签+彩色背景
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
  
  labs(title = "PCA projection based on LPS: IL-1RA, IL-1β, TNF, IL-6",
       x = paste0("PC1 (", round(summary(pca_res)$importance[2,1]*100, 2), "%)"),
       y = paste0("PC2 (", round(summary(pca_res)$importance[2,2]*100, 2), "%)"),
       color = "Cluster / Stimulus",
       shape = "Group") +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )



# ==== 选定未用于PCA的 marker ====
other_markers <- c("CCL2", "CCL3", "CCL4", "IL_10", "IL_8_1")

# ==== 筛选 LPS + 这些 marker ====
ldata_box <- ldata %>%
  filter(stimulation == "LPS", variable %in% other_markers) %>%
  filter(!is.na(level))

# ==== 加入 PCA cluster 信息（来自 pca_df） ====
cluster_df <- pca_df[, c("ID", "cluster")]
ldata_box <- merge(ldata_box, cluster_df, by = "ID")

# ==== Boxplot 绘图 ====
ggplot(ldata_box, aes(x = cluster, y = level, fill = cluster)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75)) +
  geom_jitter(shape = 21, color = "black", alpha = 0.3,
              position = position_dodge(width = 0.75), stroke = 0.3) +
  facet_wrap(~ variable, ncol = 3, scales = "free_y") +
  scale_fill_manual(values = c("1" = "#E41A1C", "2" = "#4DAF4A")) +
  theme_bw(base_size = 14) +
  labs(
    title = "Other Marker Expression by PCA Cluster (LPS)",
    x = "Cluster", y = "Log10 Marker Level"
  ) +
  theme(
    strip.text = element_text(size = 12),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    axis.title = element_text(size = 13),
    legend.position = "none"
  )



boxplot( pca_coords_HV[,1], pca_res$x[,1][kmeansres$cluster==1], pca_res$x[,1][kmeansres$cluster==2], ylim=c(-7,7))
