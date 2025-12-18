#run the code 'CTS_all_RNA_makescore_table1' first
#this is for CAP

# run the code 'CTS_all_RNA_makescore_table1' first
# this is for CAP

## ===============================================================
## 0. Packages
## ===============================================================
suppressPackageStartupMessages({
  library(randomForest)
  library(tidyverse)
  library(pheatmap)
  library(svglite)
})

dir.create("Results", showWarnings = FALSE)

## ===============================================================
## 1. Train RF model on Sepsis (NO ComBat)
## ===============================================================

common_genes <- intersect(rownames(exp_core_g), rownames(expr_logcpm))

core_mat <- exp_core_g[common_genes, , drop = FALSE]
cap_mat  <- expr_logcpm[common_genes, , drop = FALSE]

y_core <- factor(core_samples[colnames(core_mat), "CTS"])

set.seed(123)
rf_model <- randomForest(x = t(core_mat), y = y_core)

pred_noCombat <- data.frame(
  Sample = colnames(cap_mat),
  CTS_noCombat = predict(rf_model, newdata = t(cap_mat))
)

## ===============================================================
## 2. Extract ComBat-based predictions
## ===============================================================
pred_combat <- res$predictions[, c("Sample","CTS")]
names(pred_combat) <- c("Sample","CTS_ComBat")

## ===============================================================
## 3. Merge
## ===============================================================
cmp <- inner_join(pred_noCombat, pred_combat, by = "Sample")

## ===============================================================
## 4. Plot: CTS proportions (Before vs After ComBat)
## ===============================================================
cmp_long <- cmp %>%
  pivot_longer(cols = c(CTS_noCombat, CTS_ComBat),
               names_to = "Version",
               values_to = "CTS")

cts_cols <- c("1"="royalblue","2"="#B2DF8A","3"="orange")

p_prop <- ggplot(cmp_long, aes(Version, fill = CTS)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = cts_cols) +
  theme_classic(base_size = 12) +
  labs(title = "CAP CTS subtype proportions",
       y = "Proportion", x = NULL)

ggsave("Results/CTS_Proportion_BeforeAfter_ComBat.svg", p_prop, width = 5, height = 4)
ggsave("Results/CTS_Proportion_BeforeAfter_ComBat.pdf", p_prop, width = 5, height = 4)

## ===============================================================
## 5. Confusion Matrix (No-ComBat vs ComBat)
## ===============================================================
cmp <- cmp %>%
  mutate(
    CTS_noCombat = factor(CTS_noCombat, levels = c("1","2","3")),
    CTS_ComBat   = factor(CTS_ComBat,   levels = c("1","2","3"))
  )

tab_change <- table(NoComBat = cmp$CTS_noCombat,
                    ComBat   = cmp$CTS_ComBat)

write.csv(as.data.frame(tab_change),
          "Results/CTS_change_confusion_matrix.csv",
          row.names = FALSE)

## ---- Heatmap (PDF, SVG, PNG) with nicer colours ----

# high-contrast blue palette
col_heat <- colorRampPalette(c("#ffffff", "#c6dbef", "#6baed6", "#08519c"))(100)

# PDF
pheatmap(
  tab_change,
  display_numbers = TRUE,
  number_color = "black",
  fontsize_number = 12,
  color = col_heat,
  cluster_rows = FALSE, cluster_cols = FALSE,
  border_color = NA,
  main = "CTS changes: No-ComBat vs ComBat",
  filename = "Results/CTS_change_heatmap.pdf",
  width = 4, height = 4
)

# SVG
svglite("Results/CTS_change_heatmap.svg", width = 4, height = 4)
pheatmap(
  tab_change,
  display_numbers = TRUE,
  number_color = "black",
  fontsize_number = 12,
  color = col_heat,
  cluster_rows = FALSE, cluster_cols = FALSE,
  border_color = NA,
  main = "CTS changes: No-ComBat vs ComBat"
)
dev.off()

# PNG
png("Results/CTS_change_heatmap.png", width = 1600, height = 1400, res = 300)
pheatmap(
  tab_change,
  display_numbers = TRUE,
  number_color = "black",
  fontsize_number = 12,
  color = col_heat,
  cluster_rows = FALSE, cluster_cols = FALSE,
  border_color = NA,
  main = "CTS changes: No-ComBat vs ComBat"
)
dev.off()

## ===============================================================
## 6. Agreement Rate
## ===============================================================
agree_rate <- mean(cmp$CTS_noCombat == cmp$CTS_ComBat)
write.csv(data.frame(agreement = agree_rate),
          "Results/CTS_agreement_before_after_ComBat.csv",
          row.names = FALSE)

cat(sprintf("Agreement = %.1f%%\n", 100 * agree_rate))
cat("All files saved to Results/\n")
