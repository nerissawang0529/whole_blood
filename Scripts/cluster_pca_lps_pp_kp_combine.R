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

# ==== 3. 筛选 KP + 四个marker ====
selected_markers <- c("IL_1RA", "IL_1beta", "TNF", "IL_6")

ldata_sel <- ldata %>%
  filter(stimulation !="M", variable %in% selected_markers) %>%
  filter(!is.na(level))

# 宽格式
ldata_wide <- reshape2::dcast(ldata_sel, ID ~ variable*stimulation, value.var = "level")
rownames(ldata_wide) <- ldata_wide$ID
ldata_wide$ID <- NULL

ldata_wide2 <- ldata_wide[complete.cases(ldata_wide),]

# ==== 4. PCA & kmeans聚类 ====
pca_res <- prcomp(ldata_wide2, scale. = TRUE)
pca_res$x[,1] <- -pca_res$x[,1]       # Flip PC1
pca_res$rotation[,1] <- -pca_res$rotation[,1]   # Flip loadings for PC1 too
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
ldata_HV <- reshape2::melt(mdata_HV, id = c("ID", "stimulation"))
ldata_HV$level <- log10(ldata_HV$value)

ldata_HV_sel <- ldata_HV %>%
  filter(stimulation != "M", variable %in% selected_markers) %>%
  filter(!is.na(level))

ldata_HV_wide <- reshape2::dcast(ldata_HV_sel, ID ~ variable*stimulation, value.var = "level")
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
loadings_df$stim <-  sapply(strsplit(loadings_df$variable, "_"), function(parts) tail(parts, 1))
loadings_df$marker <- str_remove( rownames(loadings_df), "_[^_]+$")
loadings_df$label <- loadings_df$marker






# ==== 7. 自定义颜色 ====
marker_colors <- c(
  "IL-1RA" = "#4DAF4A",
  "IL-1beta" = "#FF6347",
  "TNF" = "#FF7F00",
  "IL-6" = "#F781BF"
)



# ==== 8. 绘图 ====
arrow_scale <- 8
label_offset <- 1.1
label_size <- 5

pca_var_explained <- pca_res$sdev^2 / sum(pca_res$sdev^2)
pc1_var <- round(pca_var_explained[1] * 100, 1)
pc2_var <- round(pca_var_explained[2] * 100, 1)


ggplot(pca_all, aes(PC1, PC2, color = cluster, shape = group)) +
  geom_point(size = 2, alpha = 0.85) +
  stat_ellipse(aes(group = cluster), linewidth = 1.2) +
  
  # 载荷箭头
  # Change color = stim to linetype = stim in geom_segment
  geom_segment(data = loadings_df,
               aes(x = 0, y = 0,
                   xend = PC1 * arrow_scale,
                   yend = PC2 * arrow_scale,
                   colour  = stim),
               arrow = arrow(length = unit(0.35, "cm")),
               linewidth = 1.5,
               inherit.aes = FALSE) +
  scale_color_manual(values=c("black","darkgray","cyan4", "red","cornflowerblue","green2"  ))+
  geom_label_repel(data = loadings_df,
                   aes(x = PC1 * arrow_scale ,
                       y = PC2 * arrow_scale ,
                       label = marker, fill = stim),
                   color = "black",
                   min.segment.length = 0,
                   lineheight=2,
                   size = 4,
                   label.padding = 0.3,
                   label.size = 1,
                   segment.color = "grey30",
                   max.overlaps = Inf,
                   box.padding = 0.3,
                   force        = 2,
                   nudge_x      = -0.9,
                   direction    = "y",
                   hjust        = 1,
                   segment.size = 0.2,
                   inherit.aes = FALSE)+
  scale_fill_manual(values=c("red","cornflowerblue","green2"))+
  theme_bw()+
  labs(x = paste0("PC1 (", pc1_var, "%)"),
       y = paste0("PC2 (", pc2_var, "%)"))

  
  
  
# ==== 9. 将 cluster 标签加入原始病人数据 ====
# 提取 PCA 病人数据的 ID 和 cluster（排除 HV）
### table 1 by PC2

hist(pca_df$PC2)
summary(pca_df$PC2)

pca_df$cluster <- ifelse(pca_df$PC2 > 0.05062 , "HiPC2", "LowPC2") #median of PC2

cluster_info <- pca_df[, c("ID", "cluster")] #if want to make the box plot for PC1 AND PC2 then only run this code

# 合并 cluster 信息到 studydata_patients
studydata_patients_clustered <- merge(studydata_patients, cluster_info, by.x = "EB_id", by.y = "ID", all.x = TRUE)

# 查看结果
table(studydata_patients_clustered$cluster, useNA = "ifany")


## table
allvars <- c("inclusion_hospital", "age_yrs",	"gender",	"ethnic_group#White",	
             "BMI",	"symptom_days",
             "CCI",	"hypertension",	"cpd",	"COPD",	"diabetes",	
             "ccd","ckd",	"mneoplasm",	"immune_sup",	"cnd",	
             "qSOFA_score",	"MEWS_score",	"CURB_score",	"PSI_new",
             "crp_1_1",	 "WBC_2_1",	"Neutrophil_unit_1",	"Lymphocyte_1_1",	
             "Creatinine_value_1",
             "Platelets_value_1",	
             "Blood_Urea_Nitrogen_value_1",	"antibiotic_seven_days", 
             "length_of_stay", 
             "hospdeath",	"mortality_d30", "mortality_d90")

catvars <- c("inclusion_hospital", "gender",	"ethnic_group#White",	
             "hypertension",	"cpd",	"COPD",	"diabetes",	
             "ccd","ckd",	"mneoplasm",	"immune_sup",	"cnd",
             "antibiotic_seven_days", 
             "hospdeath",	"mortality_d30", "mortality_d90","cluster")

nonnormal <- c("age_yrs", "BMI", "symptom_days",  "CCI", "qSOFA_score",	"MEWS_score",	"CURB_score",	"PSI_new",
               "crp_1_1",	 "WBC_2_1",	"Neutrophil_unit_1",	"Lymphocyte_1_1",	
               "Creatinine_value_1",
               "Platelets_value_1",	
               "Blood_Urea_Nitrogen_value_1",	"length_of_oxygen","length_of_stay" ,"icu_stay")

tab2     <- CreateTableOne(
  vars        = allvars, 
  data        = studydata_patients_clustered, 
  strata = "cluster",
  factorVars  = catvars,
  test        = TRUE)

print(tab2, nonnormal = nonnormal, quote = TRUE, noSpaces = TRUE, smd = T, missing = T)
# Convert the CreateTableOne object to a data frame


#the following cluster is grouped by PC1 AND PC2 group####
library(ggpubr)
library(RColorBrewer)

# 设置与 PCA 一致的颜色
cluster_colors <- c("1" = "black", "2" = "darkgray")

# Boxplot + Wilcoxon test for ALL markers
ggplot(ldata_cluster, aes(x = cluster, y = level, fill = cluster)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85, color = "black") +
  geom_jitter(width = 0.2, alpha = 0.4, size = 0.7, color = "black") +
  facet_wrap(~ variable, scales = "free_y", ncol = 4) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.y.npc = "top") +
  scale_fill_manual(values = cluster_colors) +
  labs(x = "Cluster", y = "log10 Level", title = "Marker Expression by Cluster") +
  theme_bw() +
  theme(strip.text = element_text(size = 10),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        legend.position = "none")

#####




library(pheatmap)

# === 1. Combine PC2 and selected biomarkers ===
heatmap_data <- data.frame(
  PC2 = pca_df$PC2,
  TNF_LPS = ldata_wide2$TNF_LPS,
  IL_1RA_LPS = ldata_wide2$IL_1RA_LPS,
  TNF_KP = ldata_wide2$TNF_KP,
  IL_1RA_KP = ldata_wide2$IL_1RA_KP,
  TNF_PP = ldata_wide2$TNF_PP,
  IL_1RA_PP = ldata_wide2$IL_1RA_PP
)


# === 2. Order rows (patients) by PC2 ===
heatmap_data_ordered <- heatmap_data[order(heatmap_data$PC2), ]

# === 3. Extract biomarker matrix and do Z-score normalization ===
biomarker_matrix <- as.matrix(heatmap_data_ordered[, c(
  "TNF_LPS", "IL_1RA_LPS",
  "TNF_KP", "IL_1RA_KP",
  "TNF_PP", "IL_1RA_PP"
)])
biomarker_matrix <- scale(biomarker_matrix)

# Example input
# expression_matrix: rows = markers, cols = patients
# pc2_vector: numeric vector of PC2 values (length = number of patients)

# Calculate correlations for each row against PC2
cor_results <- apply(biomarker_matrix, 2, function(x) cor(x, heatmap_data_ordered$PC2 , 
                                                          method = "pearson"))

# Convert to a named data frame
cor_df <- data.frame(
  marker = names(cor_results),
  correlation_with_PC2 = cor_results
)

cor_df <- cor_df[order(-cor_df$correlation_with_PC2),]

toplot <- t(biomarker_matrix)
toplot <- toplot[row.names(cor_df),]

# === 4. Create annotation bar for PC2 ===
annotation_col <- data.frame(PC2 = heatmap_data_ordered$PC2)
rownames(annotation_col) <- rownames(heatmap_data_ordered)

# Define gradient color for PC2 bar (blue → white → red)
pc2_colors <- colorRampPalette(c("blue", "white", "red"))(100)
annotation_colors <- list(PC2 = pc2_colors)

# === 5. Generate symmetric color breaks for heatmap ===
custom_palette <- colorRampPalette(c("blue", "grey90", "red"))(200)
max_abs <- max(abs(biomarker_matrix), na.rm = TRUE)
breaks <- seq(-max_abs, max_abs, length.out = 201)

# === 6. Draw heatmap with PC2 bar on top ===
pheatmap(
  toplot,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  color = custom_palette,
  breaks = breaks,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  main = "Z-scored Cytokine Responses (LPS, KP, PP) Ordered by PC2"
)





library(pheatmap)

# Combine PC1 and selected biomarkers
heatmap_data <- data.frame(
  PC1 = pca_df$PC1,
  TNF_LPS = ldata_wide2$TNF_LPS,
  IL_1RA_LPS = ldata_wide2$IL_1RA_LPS,
  TNF_KP = ldata_wide2$TNF_KP,
  IL_1RA_KP = ldata_wide2$IL_1RA_KP,
  TNF_PP = ldata_wide2$TNF_PP,
  IL_1RA_PP = ldata_wide2$IL_1RA_PP
)

# Sort by PC1 in descending order
heatmap_data_ordered <- heatmap_data[order(heatmap_data$PC1, decreasing = TRUE), ]

# Extract biomarker matrix and apply Z-score
biomarker_matrix <- as.matrix(heatmap_data_ordered[, c(
  "TNF_LPS", "IL_1RA_LPS", "TNF_KP", "IL_1RA_KP", "TNF_PP", "IL_1RA_PP"
)])
biomarker_matrix <- scale(biomarker_matrix)

# Create top annotation bar
annotation_col <- data.frame(PC1 = heatmap_data_ordered$PC1)
rownames(annotation_col) <- rownames(heatmap_data_ordered)

# Color gradients
pc1_colors <- colorRampPalette(c("blue", "white", "red"))(100)
annotation_colors <- list(PC1 = pc1_colors)

custom_palette <- colorRampPalette(c("blue", "grey90", "red"))(200)
max_abs <- max(abs(biomarker_matrix), na.rm = TRUE)
breaks <- seq(-max_abs, max_abs, length.out = 201)

# Draw the heatmap
pheatmap(
  t(biomarker_matrix),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  color = custom_palette,
  breaks = breaks,
  main = "Z-scored Cytokine Responses Ordered by PC1 (Descending)"
)


library(pheatmap)

# === 1. Prepare the data frame with ordering variable and markers to be plotted ===
heatmap_data <- data.frame(
  TNF_LPS = ldata_wide2$TNF_LPS,         # Ordering variable (used only for sorting and top bar)
  IL_1RA_LPS = ldata_wide2$IL_1RA_LPS,
  TNF_KP = ldata_wide2$TNF_KP,
  IL_1RA_KP = ldata_wide2$IL_1RA_KP,
  TNF_PP = ldata_wide2$TNF_PP,
  IL_1RA_PP = ldata_wide2$IL_1RA_PP
)

# === 2. Order patients by TNF_LPS ===
heatmap_data_ordered <- heatmap_data[order(heatmap_data$TNF_LPS), ]

# === 3. Extract biomarker matrix (excluding TNF_LPS) and apply Z-score normalization ===
biomarker_matrix <- as.matrix(heatmap_data_ordered[, c(
  "IL_1RA_LPS", "TNF_KP", "IL_1RA_KP", "TNF_PP", "IL_1RA_PP"
)])
biomarker_matrix <- scale(biomarker_matrix)

# === 4. Create annotation bar (top of heatmap) using raw TNF_LPS values ===
annotation_col <- data.frame(TNF_LPS = heatmap_data_ordered$TNF_LPS)
rownames(annotation_col) <- rownames(heatmap_data_ordered)

# Define a diverging color gradient: low (blue) → mid (white) → high (red)
tnf_colors <- colorRampPalette(c("blue", "white", "red"))(100)
annotation_colors <- list(TNF_LPS = tnf_colors)

# === 5. Set color palette and breaks for heatmap (centered at 0) ===
custom_palette <- colorRampPalette(c("blue", "grey90", "red"))(200)
max_abs <- max(abs(biomarker_matrix), na.rm = TRUE)
breaks <- seq(-max_abs, max_abs, length.out = 201)

# === 6. Draw heatmap ===
pheatmap(
  t(biomarker_matrix),               # Markers as rows, patients as columns
  cluster_rows = FALSE,             # Do not cluster markers
  cluster_cols = FALSE,             # Do not cluster patients
  show_colnames = FALSE,            # Hide column names (patient IDs)
  annotation_col = annotation_col,  # Add top bar showing TNF_LPS gradient
  annotation_colors = annotation_colors,
  color = custom_palette,
  breaks = breaks,
  main = "Heatmap Ordered by TNF_LPS (TNF Not Displayed)"
)



#UMAP####
install.packages("uwot")
# ==== Load required packages ====
library(uwot)
library(tidyverse)
library(ggplot2)

# ==== 1. Run UMAP and keep model ====
set.seed(123)  # For reproducibility
umap_result <- umap(
  ldata_wide2,               # Your patient-wide matrix
  scale = TRUE,              # Scale variables
  n_neighbors = 15,          # UMAP parameter: local neighborhood size
  min_dist = 0.1,            # UMAP parameter: how tightly to pack the clusters
  ret_model = TRUE           # Return model for projecting new data
)

# ==== 2. Extract UMAP embedding coordinates ====
umap_embedding <- umap_result$embedding

# ==== 3. Format patient data (221 samples) ====
umap_df <- data.frame(
  UMAP1 = umap_embedding[, 1],
  UMAP2 = umap_embedding[, 2]
)
umap_df$ID <- rownames(ldata_wide2)
umap_df$cluster <- as.factor(kmeansres$cluster)  # Use same kmeans clustering result
umap_df$group <- "Patient"

# ==== 4. Project HV samples into trained UMAP model ====
umap_HV_coords <- umap_transform(ldata_HV_wide, umap_result)

umap_HV_df <- data.frame(
  UMAP1 = umap_HV_coords[, 1],
  UMAP2 = umap_HV_coords[, 2]
)
umap_HV_df$ID <- rownames(ldata_HV_wide)
umap_HV_df$group <- "Healthy group"
umap_HV_df$cluster <- "HV"  # Mark healthy volunteers as a separate cluster

# ==== 5. Combine patient and HV for joint visualization ====
umap_all <- bind_rows(umap_df, umap_HV_df)

# ==== 6. Plot UMAP results with clusters and group labels ====
ggplot(umap_all, aes(x = UMAP1, y = UMAP2, color = cluster, shape = group)) +
  geom_point(size = 2, alpha = 0.85) +                          # Scatter points
  stat_ellipse(aes(group = cluster), linewidth = 1.2) +        # Cluster boundary ellipses
  scale_color_manual(values = c("black", "darkgray", "cyan4")) +  # Custom cluster colors
  theme_bw() +
  labs(
    title = "UMAP of Cytokine Responses",
    x = "UMAP Dimension 1",
    y = "UMAP Dimension 2"
  ) +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(size = 14, face = "bold")
  )


# ==== Step 1. 加载所需包 ====
library(destiny)
library(ggplot2)
library(dplyr)

# ==== Step 2. 准备表达矩阵 ====
# 使用你已有的 ldata_wide2，行为病人，列为已 log10 转换的 marker
expression_matrix <- as.matrix(ldata_wide2)

# ==== Step 3. 计算 Diffusion Map ====
dm <- DiffusionMap(expression_matrix)

# ==== Step 4. 提取前两个 Diffusion Component ====
diff_map_df <- as.data.frame(eigenvectors(dm)[, 1:2])
colnames(diff_map_df) <- c("DC1", "DC2")
diff_map_df$ID <- rownames(expression_matrix)

# ==== Step 5. 合并 cluster 信息 ====
# 使用你已有的 umap_df 中的 cluster 信息
diff_map_df <- left_join(diff_map_df, umap_df[, c("ID", "cluster")], by = "ID")

# ==== Step 6. 绘制 Diffusion Map ====
ggplot(diff_map_df, aes(x = DC1, y = DC2, color = cluster)) +
  geom_point(size = 2, alpha = 0.9) +
  scale_color_manual(values = c("1" = "black", "2" = "darkgray")) +
  labs(
    title = "Diffusion Map of Cytokine Responses",
    x = "Diffusion Component 1",
    y = "Diffusion Component 2"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.title = element_blank()
  )

# 添加两个 marker 到 diff_map_df
diff_map_df$TNF_LPS <- ldata_wide2[diff_map_df$ID, "TNF_LPS"]
diff_map_df$IL_1RA_LPS <- ldata_wide2[diff_map_df$ID, "IL_1RA_LPS"]

# 绘图：颜色代表 TNF，点大小代表 IL-1RA
diff_map_df$DC1 <- -diff_map_df$DC1
ggplot(diff_map_df, aes(x = DC1, y = DC2)) +
  geom_point(
    aes(color = TNF_LPS, size = IL_1RA_LPS),
    alpha = 0.9
  ) +
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = median(diff_map_df$TNF_LPS, na.rm = TRUE),
    name = "log10(TNF_LPS)"
  ) +
  scale_size(
    range = c(1, 6),
    name = "log10(IL_1RA_LPS)"
  ) +
  labs(
    title = "Diffusion Map: TNF_LPS (Color) + IL_1RA_LPS (Size)",
    x = "Diffusion Component 1",
    y = "Diffusion Component 2"
  ) +
  theme_bw()

library(patchwork)

p1 <- ggplot(diff_map_df, aes(x = DC1, y = DC2, color = TNF_LPS)) +
  geom_point(size = 2) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = median(diff_map_df$TNF_LPS)) +
  labs(title = "TNF_LPS on Diffusion Map") +
  theme_bw()

p2 <- ggplot(diff_map_df, aes(x = DC1, y = DC2, color = IL_1RA_LPS)) +
  geom_point(size = 2) +
  scale_color_gradient2(low = "green", mid = "white", high = "black", midpoint = median(diff_map_df$IL_1RA_LPS)) +
  labs(title = "IL_1RA_LPS on Diffusion Map") +
  theme_bw()

# 拼图
p1 + p2
library(ggplot2)
library(patchwork)

# 先确保 diff_map_df 包含所有变量
diff_map_df$TNF_PP     <- ldata_wide2[diff_map_df$ID, "TNF_PP"]
diff_map_df$IL_1RA_PP  <- ldata_wide2[diff_map_df$ID, "IL_1RA_PP"]
diff_map_df$TNF_KP     <- ldata_wide2[diff_map_df$ID, "TNF_KP"]
diff_map_df$IL_1RA_KP  <- ldata_wide2[diff_map_df$ID, "IL_1RA_KP"]

# 如果尚未翻转 DC1，可手动执行：
diff_map_df$DC1 <- -diff_map_df$DC1

# LPS
p1 <- ggplot(diff_map_df, aes(x = DC1, y = DC2, color = TNF_LPS)) +
  geom_point(size = 2) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red",
                        midpoint = median(diff_map_df$TNF_LPS, na.rm = TRUE)) +
  labs(title = "TNF_LPS on Diffusion Map") +
  theme_bw()

p2 <- ggplot(diff_map_df, aes(x = DC1, y = DC2, color = IL_1RA_LPS)) +
  geom_point(size = 2) +
  scale_color_gradient2(low = "green", mid = "white", high = "black",
                        midpoint = median(diff_map_df$IL_1RA_LPS, na.rm = TRUE)) +
  labs(title = "IL_1RA_LPS on Diffusion Map") +
  theme_bw()

# PP
p3 <- ggplot(diff_map_df, aes(x = DC1, y = DC2, color = TNF_PP)) +
  geom_point(size = 2) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red",
                        midpoint = median(diff_map_df$TNF_PP, na.rm = TRUE)) +
  labs(title = "TNF_PP on Diffusion Map") +
  theme_bw()

p4 <- ggplot(diff_map_df, aes(x = DC1, y = DC2, color = IL_1RA_PP)) +
  geom_point(size = 2) +
  scale_color_gradient2(low = "green", mid = "white", high = "black",
                        midpoint = median(diff_map_df$IL_1RA_PP, na.rm = TRUE)) +
  labs(title = "IL_1RA_PP on Diffusion Map") +
  theme_bw()

# KP
p5 <- ggplot(diff_map_df, aes(x = DC1, y = DC2, color = TNF_KP)) +
  geom_point(size = 2) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red",
                        midpoint = median(diff_map_df$TNF_KP, na.rm = TRUE)) +
  labs(title = "TNF_KP on Diffusion Map") +
  theme_bw()

p6 <- ggplot(diff_map_df, aes(x = DC1, y = DC2, color = IL_1RA_KP)) +
  geom_point(size = 2) +
  scale_color_gradient2(low = "green", mid = "white", high = "black",
                        midpoint = median(diff_map_df$IL_1RA_KP, na.rm = TRUE)) +
  labs(title = "IL_1RA_KP on Diffusion Map") +
  theme_bw()

# 拼图输出为 3 行 2 列
(p1 | p2) / (p3 | p4) / (p5 | p6)


# 将 DC1 加入为 pseudotime
pseudotime_df <- diff_map_df[, c("ID", "DC1")]
colnames(pseudotime_df)[2] <- "pseudotime"
studydata_patients_pseudo <- merge(studydata_patients, pseudotime_df, by.x = "EB_id", by.y = "ID", all.x = TRUE)

# 三分组 cut
studydata_patients_pseudo$pseudo_group <- cut(
  studydata_patients_pseudo$pseudotime,
  breaks = quantile(studydata_patients_pseudo$pseudotime, probs = c(0, 0.33, 0.66, 1), na.rm = TRUE),
  labels = c("Low", "Mid", "High")
)

# 示例：死亡率比较
table(studydata_patients_pseudo$pseudo_group, studydata_patients_pseudo$mortality_d30)

ggplot(diff_map_df, aes(x = DC1, y = TNF_LPS)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE, color = "red")

