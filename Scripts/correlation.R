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

# === Load required libraries ===
library(tidyverse)
library(Hmisc)
library(ggcorrplot)

# === Compute correlation and p-values ===
cor_res <- rcorr(as.matrix(ldata_wide2))
cor_matrix <- cor_res$r
p_matrix <- cor_res$P

# === Set non-significant correlations to NA for coloring ===
sig_level <- 0.05
masked_cor_matrix <- cor_matrix
masked_cor_matrix[p_matrix >= sig_level] <- NA  # nonsignificant values become NA

library(tidyverse)
library(Hmisc)
library(ggcorrplot)

# === Step 1: Compute correlation and p-values ===
cor_res <- rcorr(as.matrix(ldata_wide2))
cor_matrix <- cor_res$r
p_matrix <- cor_res$P

# === Step 2: Reorder matrix according to desired marker order ===
desired_order <- c(
  "TNF_LPS", "TNF_PP", "TNF_KP",
  "IL_1RA_LPS", "IL_1RA_PP", "IL_1RA_KP",
  "IL_1beta_LPS", "IL_1beta_PP", "IL_1beta_KP",
  "IL_6_LPS", "IL_6_KP", "IL_6_PP"
)

# Ensure the order exists in your data
desired_order <- intersect(desired_order, colnames(cor_matrix))

# Reorder correlation and p-value matrices
cor_matrix <- cor_matrix[desired_order, desired_order]
p_matrix <- p_matrix[desired_order, desired_order]

# === Step 3: Plot with fixed color scale from 0 to 1 ===
ggcorrplot(
  cor_matrix,
  method = "circle",
  type = "lower",
  lab = TRUE,
  lab_size = 3,
  outline.color = "white",
  colors = c("white", "red"),        # Only positive correlation
  show.legend = TRUE,
  legend.title = "Corr",
  tl.cex = 10,
  p.mat = p_matrix,
  sig.level = 0.05,
  insig = "blank",
  hc.order = FALSE,
  ggtheme = theme_minimal()
) +
  scale_fill_gradient(
    limits = c(0, 1),
    low = "white", high = "red",
    name = "Corr",
    breaks = c(0, 0.5, 1)
  ) +
  ggtitle("Correlation Heatmap") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 16)
  )
