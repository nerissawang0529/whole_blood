#pca cluster:the original data

rm(list = ls())

studydata_patients <- read.csv("Original_data/studydata_patients.csv")
#studydata_HV <- read.csv("Original_data/studydata_HV.csv")
#studydata_combined <- rbind(studydata_patients, studydata_HV)

Luminex_CAP_COVID_HV <- read.csv("Original_data/Luminex_CAP_COVID_HV.csv")
merged_data <- merge(Luminex_CAP_COVID_HV, studydata_patients, 
                     by.x = "ID", by.y = "EB_id")



mdata <- merged_data[,1:12]
mdata$Day <- NULL

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


resid_w["3005", "PP_IL_1beta"] <- 0
resid_w["3350", "PP_CCL3"] <- 0


pca_res <- prcomp(resid_w)

library(ggfortify)

autoplot(pca_res)


autoplot(pca_res, 
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3)


kmeansres <- kmeans(pca_res$x[,1:2], centers = 2) 

pat_clust <- data.frame(cluster=as.factor(kmeansres$cluster))

autoplot(pca_res, 
         data = pat_clust, 
         colour = "cluster",
         loadings = TRUE, 
         loadings.colour = 'blue',
         loadings.label = TRUE, 
         loadings.label.size = 3,
         frame = TRUE, 
         frame.type = 'norm')  # or 't' or 'euclid' depending on the look you want


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
  facet_wrap(~ stimulation * variable, ncol = 9, scales = "free_y") +  
  scale_fill_manual(values = c("1" = "#E41A1C", "2" = "#4DAF4A")) +  # PCA颜色同步
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
    title = "Profile of markers"
  )




# Load and process HV data
studydata_HV <- read.csv("Original_data/studydata_HV.csv")
Luminex_HV <- read.csv("Original_data/Luminex_CAP_COVID_HV.csv")

merged_HV <- merge(Luminex_HV, studydata_HV, by.x = "ID", by.y = "EB_id")
mdata_HV <- merged_HV[,1:12]
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

# ====== PCA loading 数据准备 ======
loadings_df <- as.data.frame(pca_res$rotation[, 1:2])
loadings_df$variable <- rownames(loadings_df)
loadings_df$stim <- gsub("_.*", "", loadings_df$variable)
loadings_df$marker <- gsub("^[^_]+_", "", loadings_df$variable)
loadings_df$marker <- gsub("_", "-", loadings_df$marker)
loadings_df$label <- loadings_df$marker

# ====== marker颜色，突出 IL-8-1 和 IL-1beta ======
unique_markers <- unique(loadings_df$marker)
marker_colors <- setNames(
  c(
    "#E41A1C",  # CCL2
    "#FFFF33",  # CCL3
    "#377EB8",  # CCL4
    "#666666",  # IL-10
    "#FF6347",  # IL-1beta 🔸 NEW COLOR
    "#4DAF4A",  # IL-1RA
    "#F781BF",  # IL-6
    "#6A5ACD",  # IL-8-1 🔸 NEW COLOR
    "#FF7F00"   # TNF
  ),
  c("CCL2", "CCL3", "CCL4", "IL-10", "IL-1beta", "IL-1RA", "IL-6", "IL-8-1", "TNF")
)

# 更新 group 显示
pca_all$group <- ifelse(pca_all$group == "HV", "Healthy group", pca_all$group)

# 箭头颜色（刺激来源）
stim_cols_arrow <- c("KP" = "#984EA3", "LPS" = "#FF7F00", "PP" = "#A65628")

# 参数控制
arrow_scale <- 16       # 🔺更长箭头
label_offset <- 1.1
label_size <- 5.2

# ====== 绘图主体 ======
ggplot(pca_all, aes(PC1, PC2, color = cluster, shape = group)) +
  geom_point(size = 2, alpha = 0.85) +
  stat_ellipse(aes(group = cluster), linewidth = 1.2) +
  
  # 箭头
  geom_segment(data = loadings_df,
               aes(x = 0, y = 0,
                   xend = PC1 * arrow_scale,
                   yend = PC2 * arrow_scale,
                   color = stim),
               arrow = arrow(length = unit(0.35, "cm")),
               linewidth = 1.2,
               inherit.aes = FALSE) +
  
  # 标签 + 彩色背景（marker）
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
  
  # 自定义颜色
  scale_color_manual(
    values = c("1" = "#E41A1C", "2" = "#4DAF4A", "HV" = "#377EB8", stim_cols_arrow),
    breaks = c("1", "2", "HV", "KP", "LPS", "PP"),
    labels = c("Cluster 1", "Cluster 2", "Healthy group", "kneu stimulus", "LPS stimulus", "pneu stimulus")
  ) +
  scale_fill_manual(values = marker_colors, name = "Marker") +
  scale_shape_manual(values = c("Patient" = 17, "Healthy group" = 16)) +
  
  # 标签与主题
  labs(title = "PCA projection with All Cytokine Loadings",
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

ggplot(pca_all[pca_all$group == "Patient", ], 
       aes(PC1, PC2, color = cluster)) +
  geom_point(size = 2, alpha = 0.85) +
  stat_ellipse(aes(group = cluster), linewidth = 1.2) +
  scale_color_manual(
    values = c("1" = "#E41A1C", "2" = "#4DAF4A"),
    labels = c("Cluster 1", "Cluster 2")
  ) +
  theme_minimal(base_size = 14) +
  labs(title = "PCA projection of Patient Clusters",
       x = paste0("PC1 (", round(summary(pca_res)$importance[2,1]*100, 2), "%)"),
       y = paste0("PC2 (", round(summary(pca_res)$importance[2,2]*100, 2), "%)"),
       color = "Cluster")

#add HV project prodict from pca prcomp in r
#residuals to cluster is ok
#try raw data
#no medium use HV as ctr; media of HV/CAP in the box plot lines

#step 1
ggplot(pca_all[pca_all$group == "Patient", ], 
       aes(PC1, PC2, color = cluster)) +
  geom_point(size = 2, alpha = 0.85) +
  stat_ellipse(aes(group = cluster), linewidth = 1.2) +
  scale_color_manual(
    values = c("1" = "#E41A1C", "2" = "#4DAF4A"),
    labels = c("Cluster 1", "Cluster 2")
  ) +
  theme_minimal(base_size = 14) +
  labs(title = "PCA projection of Patient Clusters",
       x = paste0("PC1 (", round(summary(pca_res)$importance[2,1]*100, 2), "%)"),
       y = paste0("PC2 (", round(summary(pca_res)$importance[2,2]*100, 2), "%)"),
       color = "Cluster")

#step 2
ggplot(pca_all, 
       aes(PC1, PC2, color = cluster, shape = group)) +
  geom_point(size = 2, alpha = 0.85) +
  stat_ellipse(aes(group = cluster), linewidth = 1.2) +
  scale_color_manual(
    values = c("1" = "#E41A1C", "2" = "#4DAF4A", "HV" = "#377EB8"),
    labels = c("Cluster 1", "Cluster 2", "Healthy group")
  ) +
  scale_shape_manual(values = c("Patient" = 17, "Healthy group" = 16)) +
  theme_minimal(base_size = 14) +
  labs(title = "PCA: Patients and Healthy Volunteers",
       x = paste0("PC1 (", round(summary(pca_res)$importance[2,1]*100, 2), "%)"),
       y = paste0("PC2 (", round(summary(pca_res)$importance[2,2]*100, 2), "%)"),
       color = "Cluster",
       shape = "Group")

#step 3
# 缩短箭头长度
arrow_scale <- 4  # ⬅️ 再次缩短箭头

ggplot(pca_all, aes(PC1, PC2, color = cluster, shape = group)) +
  geom_point(size = 2, alpha = 0.85) +
  stat_ellipse(aes(group = cluster), linewidth = 1.2) +
  
  # 箭头（刺激来源）
  geom_segment(data = loadings_df,
               aes(x = 0, y = 0,
                   xend = PC1 * arrow_scale,
                   yend = PC2 * arrow_scale,
                   color = stim),
               arrow = arrow(length = unit(0.3, "cm")),
               linewidth = 1.2,
               inherit.aes = FALSE) +
  
  scale_color_manual(
    values = c("1" = "#E41A1C", "2" = "#4DAF4A", "HV" = "#377EB8", 
               "KP" = "#984EA3", "LPS" = "#FF7F00", "PP" = "#A65628"),
    labels = c("Cluster 1", "Cluster 2", "Healthy group", 
               "kneu stimulus", "LPS stimulus", "pneu stimulus")
  ) +
  scale_shape_manual(values = c("Patient" = 17, "Healthy group" = 16)) +
  
  labs(title = "PCA with Cytokine Arrows (Stimulus Only)",
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


#step 4
ggplot(pca_all, aes(PC1, PC2, color = cluster, shape = group)) +
  geom_point(size = 2, alpha = 0.85) +
  stat_ellipse(aes(group = cluster), linewidth = 1.2) +
  
  geom_segment(data = loadings_df,
               aes(x = 0, y = 0,
                   xend = PC1 * 16,
                   yend = PC2 * 16,
                   color = stim),
               arrow = arrow(length = unit(0.35, "cm")),
               linewidth = 1.2,
               inherit.aes = FALSE) +
  
  geom_label_repel(data = loadings_df,
                   aes(x = PC1 * 16 * 1.1,
                       y = PC2 * 16 * 1.1,
                       label = label, fill = marker),
                   color = "black",
                   size = 4.2,
                   label.padding = 0.3,
                   label.size = 0.4,
                   segment.color = "grey30",
                   max.overlaps = Inf,
                   box.padding = 0.5,
                   inherit.aes = FALSE) +
  
  scale_color_manual(
    values = c("1" = "#E41A1C", "2" = "#4DAF4A", "HV" = "#377EB8", 
               "KP" = "#984EA3", "LPS" = "#FF7F00", "PP" = "#A65628"),
    labels = c("Cluster 1", "Cluster 2", "Healthy group", 
               "kneu stimulus", "LPS stimulus", "pneu stimulus")
  ) +
  scale_fill_manual(values = marker_colors, name = "Marker") +
  scale_shape_manual(values = c("Patient" = 17, "Healthy group" = 16)) +
  theme_minimal(base_size = 14) +
  labs(title = "PCA projection with All Cytokine Loadings",
       x = paste0("PC1 (", round(summary(pca_res)$importance[2,1]*100, 2), "%)"),
       y = paste0("PC2 (", round(summary(pca_res)$importance[2,2]*100, 2), "%)"),
       color = "Cluster / Stimulus",
       shape = "Group")

