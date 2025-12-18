# 安装一次即可
install.packages("tableone")

# 每次使用前加载
library(tableone)


#pca cluster:the original data

rm(list = ls())

studydata_patients <- read.csv("Original_data/studydata_patients.csv")
#studydata_HV <- read.csv("Original_data/studydata_HV.csv")
#studydata_combined <- rbind(studydata_patients, studydata_HV)

Luminex_CAP_COVID_HV <- read.csv("Original_data/Luminex_CAP_COVID_HV.csv")
merged_data <- merge(Luminex_CAP_COVID_HV, studydata_patients, 
                     by.x = "ID", by.y = "EB_id")



mdata <- merged_data[, c(1, 12,14:22)]
mdata <- mdata[!apply(is.na(mdata[, -c(1, 2)]), 1, all), ]
mdata[, -1] <- lapply(mdata[, -1], function(x) ifelse(x < 0, 1, x))


library(reshape2)
ldata <- reshape2::melt(mdata, id=c("ID","stimulation") )
ldata$level <- log10(ldata$value)
ldata <- ldata[ldata$stimulation != "M",]

hist(ldata$level)
library(ggplot2)
ggplot(ldata, aes(x=level, fill=stimulation))+
  geom_histogram(alpha=0.8, binwidth = 0.2)+
  facet_wrap(~variable)+
  theme_classic()

library(lme4)
library(lmerTest)

ldata$stimulation <- as.factor(ldata$stimulation)
ldata$stimulation <- factor(ldata$stimulation, levels = c("KP","LPS","PP"))

ldata_keep <- ldata

ldata <- ldata[!is.na(ldata$level),]


mod0 <- lmerTest::lmer(level ~ variable + stimulation + (1|ID) , data=ldata)

summary(mod0)


mod <- lmerTest::lmer(level ~ variable * stimulation + (1|ID) , data=ldata)

summary(mod)


anova(mod0, mod)


residualdf <- residuals(mod)

predictdf <- predict(mod)


predictdf_fixedonly <- predict(mod, re.form=NA, newdata=ldata_keep )
residualdf_fixedonly <- data.frame(ldata_keep, predicted=predictdf_fixedonly)

residualdf_fixedonly[residualdf_fixedonly$ID==3005 & 
                       residualdf_fixedonly$variable=="IL_1beta" &
                       residualdf_fixedonly$stimulation=="PP",]

residualdf_fixedonly[residualdf_fixedonly$ID==3005 & 
                       residualdf_fixedonly$variable=="IL_1beta" &
                       residualdf_fixedonly$stimulation=="PP",]

residualdf_fixedonly[is.na(residualdf_fixedonly$value),]



residualdf <- data.frame(ldata, predicted=predict(mod),  resid= residuals(mod))

ggplot(residualdf, aes(x=resid, fill=stimulation))+
  geom_histogram(alpha=0.8)+
  facet_wrap(~variable)+
  theme_classic()


#### make wide residuals 

residualdf0 <- residualdf
residualdf0$level <- NULL
residualdf0$predicted <- NULL


resid_w <- reshape2::dcast(residualdf0, ID ~ stimulation + variable, value.var="resid" )

row.names(resid_w) <- resid_w$ID
resid_w$ID <- NULL
resid_w <- resid_w[complete.cases(resid_w) & apply(resid_w, 1, function(x) all(is.finite(x))), ]

pca_res <- prcomp(resid_w)

library(ggfortify)

autoplot(pca_res)


autoplot(pca_res, 
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3)


kmeansres <- kmeans(pca_res$x[,1:3], centers = 3) 

pat_clust <- data.frame(cluster=as.factor(kmeansres$cluster))

autoplot(pca_res, data=pat_clust, colour="cluster",
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3)

autoplot(pca_res, data=pat_clust, colour="cluster")
table(pat_clust$cluster)
pat_clust$ID <- rownames(pat_clust)

##table1
studydata_patients <- read.csv("Original_data/studydata_patients.csv")
studydata_patients_0 <- merge(studydata_patients,pat_clust,by.x = "EB_id", by.y = "ID")

## table
allvars <- c("inclusion_hospital", "sampling_time",		"age_yrs",	"gender",	"ethnic_group#White",	
             "BMI",	"symptom_days",
             "CCI",	"hypertension",	"cpd",	"COPD",	"diabetes",	
             "ccd","ckd",	"mneoplasm",	"immune_sup",	"cnd",	
             "qSOFA_score",	"MEWS_score",	"CURB_score",	"PSI_new",
             "crp_1_1",	 "WBC_2_1",	"Neutrophil_unit_1",	"Lymphocyte_1_1",	
             "Creatinine_value_1",
             "Platelets_value_1",	
             "Blood_Urea_Nitrogen_value_1",	
             "Oxygen_therapy_1", 
             "length_of_oxygen", "antibiotic_seven_days", 
             "length_of_stay", 
             "ICU_Medium_Care_admission_1",	 "icu_stay", 
             "Non_invasive_ventilation_1",	"Invasive_ventilation_1",	"lenght_of_intubation",
             "bacteremia", "pathogen_cultured",
             "hospdeath",	"mortality_d30", "mortality_d90","pathogen","Gram_negative","Gram_positive","Virus","Fungi","cluster")

catvars <- c("inclusion_hospital", "gender",	"ethnic_group#White",	
             "hypertension",	"cpd",	"COPD",	"diabetes",	
             "ccd","ckd",	"mneoplasm",	"immune_sup",	"cnd",	
             "Oxygen_therapy_1", 
             "antibiotic_seven_days", 
             "ICU_Medium_Care_admission_1",	 
             "Non_invasive_ventilation_1",	"Invasive_ventilation_1",	
             "bacteremia",
             "hospdeath",	"mortality_d30", "mortality_d90", "pathogen_cultured","pathogen","Gram_negative","Gram_positive","Virus","Fungi")

nonnormal <- c("sampling_time",	"age_yrs", "BMI", "symptom_days",  "CCI", "qSOFA_score",	"MEWS_score",	"CURB_score",	"PSI_new",
               "crp_1_1",	 "WBC_2_1",	"Neutrophil_unit_1",	"Lymphocyte_1_1",	
               "Creatinine_value_1",
               "Platelets_value_1",	
               "Blood_Urea_Nitrogen_value_1",	"length_of_oxygen","length_of_stay" ,"icu_stay", "lenght_of_intubation")

tab2     <- CreateTableOne(
  vars        = allvars, 
  data        = studydata_patients_0, 
  strata = "cluster",
  factorVars  = catvars,
  test        = TRUE)

print(tab2, nonnormal = nonnormal, quote = TRUE, noSpaces = TRUE, smd = T, missing = T)




#box plot for clusters (3,9 box)
ldatamerge <- merge(ldata,pat_clust, by = "ID")
ggplot(ldatamerge, aes(x = cluster, y = level, fill = cluster)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75)) + 
  geom_jitter(shape = 21, color = "black", position = position_dodge(width = 0.75), 
              alpha = 0.3, stroke = 0.3) + 
  facet_wrap(~ stimulation*variable, ncol = 9, scales = "free_y") +  
  scale_fill_brewer(palette = "Set2") +  
  theme_bw(base_size = 16) +  
  theme(
    strip.text = element_text(size = 10),  
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12), 
    axis.text.y = element_text(size = 12), 
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 12)
  ) +
  labs(
    x = "Cluster",
    y = "Marker Value (log10 transformed)",
    title = "Profile of Non-infected controls"
  )


# Load and process HV data
studydata_HV <- read.csv("Original_data/studydata_HV.csv")
Luminex_HV <- read.csv("Original_data/Luminex_CAP_COVID_HV.csv")

merged_HV <- merge(Luminex_HV, studydata_HV, by.x = "ID", by.y = "EB_id")
mdata_HV <- merged_HV[, c(1, 12,14:22)]
mdata_HV$Day <- NULL

ldata_HV <- reshape2::melt(mdata_HV, id=c("ID","stimulation") )
ldata_HV$level <- log10(ldata_HV$value)
ldata_HV <- ldata_HV[ldata_HV$stimulation != "M",]

# Predict fixed effects for HV using patient model
ldata_HV$stimulation <- factor(ldata_HV$stimulation, levels = c("KP","LPS","PP"))
pred_HV <- predict(mod, newdata = ldata_HV, re.form = NA)
resid_HV <- ldata_HV$level - pred_HV
ldata_HV$resid <- resid_HV


resid_HV_wide <- reshape2::dcast(ldata_HV, ID ~ stimulation + variable, value.var = "resid")

# Make sure column order matches `resid_w`
resid_HV_wide <- resid_HV_wide[, colnames(resid_w)]

# Use same PCA rotation from patient-only PCA
pca_coords_HV <- predict(pca_res, newdata = resid_HV_wide)

resid_HV_wide <- reshape2::dcast(ldata_HV, ID ~ stimulation + variable, value.var = "resid")

# Make sure column order matches `resid_w`
resid_HV_wide <- resid_HV_wide[, colnames(resid_w)]

# Use same PCA rotation from patient-only PCA
pca_coords_HV <- predict(pca_res, newdata = resid_HV_wide)

pca_pat_df <- data.frame(pca_res$x[,1:2])
pca_pat_df$group <- "Patient"
pca_pat_df$ID <- rownames(pca_pat_df)

pca_HV_df <- data.frame(pca_coords_HV[,1:2])
pca_HV_df$group <- "HV"
pca_HV_df$ID <- rownames(pca_HV_df)

pca_all <- rbind(pca_pat_df, pca_HV_df)

pca_all <- merge(pca_all, pat_clust, by = "ID", all.x = TRUE)  # Cluster only exists for patients
pca_all$cluster <- ifelse(is.na(pca_all$cluster), "HV", as.character(pca_all$cluster))

library(ggplot2)
ggplot(pca_all, aes(x = PC1, y = PC2, color = cluster, shape = group)) +
  geom_point(size = 2, alpha = 0.8) +
  stat_ellipse(aes(group = cluster), linetype = "dashed", level = 0.8) +
  theme_minimal(base_size = 14) +
  labs(title = "PCA projection of Patients and Healthy Volunteers",
       x = paste0("PC1 (", round(summary(pca_res)$importance[2,1]*100, 2), "%)"),
       y = paste0("PC2 (", round(summary(pca_res)$importance[2,2]*100, 2), "%)"))

library(ggplot2)
library(ggrepel)
library(dplyr)
library(RColorBrewer)

# ==== 准备 loadings_df 数据 ====
loadings_df <- as.data.frame(pca_res$rotation[, 1:2])
loadings_df$variable <- rownames(loadings_df)
loadings_df$stim <- gsub("_.*", "", loadings_df$variable)
loadings_df$marker <- gsub("^[^_]+_", "", loadings_df$variable)
loadings_df$marker <- gsub("_", "-", loadings_df$marker)
loadings_df$label <- loadings_df$variable  # or use just marker name: loadings_df$marker

# ==== marker 颜色区分（强对比色）====
unique_markers <- unique(loadings_df$marker)
marker_palette <- RColorBrewer::brewer.pal(n = 9, "Set1")
if (length(unique_markers) > length(marker_palette)) {
  marker_palette <- colorRampPalette(marker_palette)(length(unique_markers))
}
marker_colors <- setNames(marker_palette, unique_markers)

# ==== 箭头颜色 ====
stim_cols_arrow <- c("KP" = "#984EA3", "LPS" = "#FF7F00", "PP" = "#A65628")

# ==== PCA 点位信息准备 ====
pca_all$group <- ifelse(pca_all$group == "HV", "Healthy group", pca_all$group)

# ==== 参数 ====
arrow_scale <- 12
label_offset <- 1.1
label_size <- 4     # ✅ 缩小了标签字体
label_box_padding <- 0.25
label_text_color <- "black"

# ==== 绘图 ====
ggplot(pca_all, aes(PC1, PC2, color = cluster, shape = group)) +
  geom_point(size = 2, alpha = 0.85) +
  stat_ellipse(aes(group = cluster), linewidth = 1.1) +
  
  # 箭头
  geom_segment(data = loadings_df,
               aes(x = 0, y = 0,
                   xend = PC1 * arrow_scale,
                   yend = PC2 * arrow_scale,
                   color = stim),
               arrow = arrow(length = unit(0.35, "cm")),
               linewidth = 1.1,
               inherit.aes = FALSE) +
  
  # 标签
  geom_label_repel(data = loadings_df,
                   aes(x = PC1 * arrow_scale * label_offset,
                       y = PC2 * arrow_scale * label_offset,
                       label = marker, fill = marker),
                   color = label_text_color,
                   size = label_size,
                   label.padding = label_box_padding,
                   label.size = 0.3,
                   segment.color = "grey40",
                   max.overlaps = Inf,
                   box.padding = 0.45,
                   inherit.aes = FALSE) +
  
  # 调色板
  scale_color_manual(
    values = c("1" = "#E41A1C", "2" = "#4DAF4A", "HV" = "#377EB8", stim_cols_arrow),
    breaks = c("1", "2", "HV", "KP", "LPS", "PP"),
    labels = c("Cluster 1", "Cluster 2", "Healthy group", "kneu stimulus", "LPS stimulus", "pneu stimulus")
  ) +
  scale_fill_manual(values = marker_colors, name = "Marker") +
  scale_shape_manual(values = c("Patient" = 17, "Healthy group" = 16)) +
  
  # 标签与格式
  labs(title = "PCA projection with All Cytokine Loadings",
       x = paste0("PC1 (", round(summary(pca_res)$importance[2, 1] * 100, 2), "%)"),
       y = paste0("PC2 (", round(summary(pca_res)$importance[2, 2] * 100, 2), "%)"),
       color = "Cluster / Stimulus",
       shape = "Group") +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )





#step1
ggplot(pca_all[pca_all$group == "Patient", ], 
       aes(PC1, PC2, color = cluster, shape = group)) +
  geom_point(size = 2, alpha = 0.85) +
  stat_ellipse(aes(group = cluster), linewidth = 1.1) +
  scale_color_manual(values = c("1" = "#E41A1C",  # red
                                "2" = "#4DAF4A",  # green
                                "3" = "#999999")) +  # gray
  scale_shape_manual(values = c("Patient" = 17)) +
  labs(title = "Step 1: PCA - Patient Clusters Only") +
  theme_minimal(base_size = 14)

#step2
ggplot(pca_all, aes(PC1, PC2, color = cluster, shape = group)) +
  geom_point(size = 2, alpha = 0.85) +
  stat_ellipse(aes(group = cluster), linewidth = 1.1) +
  scale_color_manual(values = c("1" = "#E41A1C", 
                                "2" = "#4DAF4A", 
                                "3" = "#999999", 
                                "HV" = "#377EB8")) +
  scale_shape_manual(values = c("Patient" = 17, "Healthy group" = 16)) +
  labs(title = "Step 2: PCA - Add Healthy Volunteers") +
  theme_minimal(base_size = 14)

#step3
ggplot(pca_all, aes(PC1, PC2, color = cluster, shape = group)) +
  geom_point(size = 2, alpha = 0.85) +
  stat_ellipse(aes(group = cluster), linewidth = 1.1) +
  geom_segment(data = loadings_df,
               aes(x = 0, y = 0,
                   xend = PC1 * 6,  # shorter arrow
                   yend = PC2 * 6,
                   color = stim),
               arrow = arrow(length = unit(0.3, "cm")),
               linewidth = 1.1,
               inherit.aes = FALSE) +
  scale_color_manual(values = c("1" = "#E41A1C", 
                                "2" = "#4DAF4A", 
                                "3" = "#999999", 
                                "HV" = "#377EB8", 
                                "KP" = "#984EA3", 
                                "LPS" = "#FF7F00", 
                                "PP" = "#A65628")) +
  scale_shape_manual(values = c("Patient" = 17, "Healthy group" = 16)) +
  labs(title = "Step 3: PCA with Cytokine Arrows (Stimulus Only)") +
  theme_minimal(base_size = 14)

#step4
ggplot(pca_all, aes(PC1, PC2, color = cluster, shape = group)) +
  geom_point(size = 2, alpha = 0.85) +
  stat_ellipse(aes(group = cluster), linewidth = 1.1) +
  
  # Arrows
  geom_segment(data = loadings_df,
               aes(x = 0, y = 0,
                   xend = PC1 * 12,
                   yend = PC2 * 12,
                   color = stim),
               arrow = arrow(length = unit(0.35, "cm")),
               linewidth = 1.1,
               inherit.aes = FALSE) +
  
  # Labels
  geom_label_repel(data = loadings_df,
                   aes(x = PC1 * 12 * 1.1,
                       y = PC2 * 12 * 1.1,
                       label = marker, fill = marker),
                   color = "black",
                   size = 4,
                   label.padding = 0.25,
                   label.size = 0.3,
                   segment.color = "grey40",
                   max.overlaps = Inf,
                   box.padding = 0.45,
                   inherit.aes = FALSE) +
  
  scale_color_manual(values = c("1" = "#E41A1C", 
                                "2" = "#4DAF4A", 
                                "3" = "#999999", 
                                "HV" = "#377EB8", 
                                "KP" = "#984EA3", 
                                "LPS" = "#FF7F00", 
                                "PP" = "#A65628")) +
  scale_fill_manual(values = marker_colors, name = "Marker") +
  scale_shape_manual(values = c("Patient" = 17, "Healthy group" = 16)) +
  
  labs(title = "Step 4: PCA with All Cytokine Loadings",
       x = paste0("PC1 (", round(summary(pca_res)$importance[2, 1] * 100, 2), "%)"),
       y = paste0("PC2 (", round(summary(pca_res)$importance[2, 2] * 100, 2), "%)"),
       color = "Cluster / Stimulus",
       shape = "Group") +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
