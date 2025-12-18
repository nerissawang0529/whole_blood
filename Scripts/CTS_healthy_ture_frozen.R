#run the code 'healthy_ran_CTS'

## ===============================================================
##  FROZEN TRANSFORM for Healthy cohort (non-infectious)
## ===============================================================

# 1. Common genes
common_genes <- intersect(rownames(exp_core_g), rownames(expr_logcpm_hv))

core_mat <- exp_core_g[common_genes, , drop = FALSE]   # sepsis reference (already ComBat)
hv_mat   <- expr_logcpm_hv[common_genes, , drop = FALSE]  # healthy raw logCPM

# 2. Gene-wise frozen mean/SD from sepsis
mu_core <- rowMeans(core_mat)
sd_core <- apply(core_mat, 1, sd)
sd_core[sd_core == 0 | is.na(sd_core)] <- 1

# 3. Gene-wise batch statistics from healthy
mu_hv <- rowMeans(hv_mat)
sd_hv <- apply(hv_mat, 1, sd)
sd_hv[sd_hv == 0 | is.na(sd_hv)] <- 1

# 4. FROZEN TRANSFORM: map healthy → frozen sepsis space
hv_frozen <- sweep(hv_mat, 1, mu_hv, "-")        # subtract HV mean
hv_frozen <- sweep(hv_frozen, 1, sd_hv, "/")     # divide by HV sd
hv_frozen <- sweep(hv_frozen, 1, sd_core, "*")   # multiply sepsis sd
hv_frozen <- sweep(hv_frozen, 1, mu_core, "+")   # add sepsis mean

# 5. Train classifier on sepsis reference (NO ComBat on healthy)
y_core <- factor(core_samples[colnames(core_mat), "CTS"])

set.seed(123)
rf_model_frozen <- randomForest(x = t(core_mat), y = y_core)

# 6. Predict CTS on frozen Healthy matrix
pred_hv_frozen <- data.frame(
  Sample = colnames(hv_frozen),
  CTS_frozen = predict(rf_model_frozen, newdata = t(hv_frozen))
)

# 7. Save output
write.csv(pred_hv_frozen,
          "Figure/cts_yourcohort_healthy_predictions_frozen.csv",
          row.names = FALSE)

cat("✅ Saved FROZEN-transform Healthy CTS to: Figure/cts_yourcohort_healthy_predictions_frozen.csv\n")


# Load frozen results
frozen <- read.csv("Figure/cts_yourcohort_healthy_predictions_frozen.csv")
frozen$CTS_frozen <- factor(frozen$CTS_frozen, levels = c("1","2","3"))

# merge 3 versions
df3 <- df %>%     # df already has before + after
  inner_join(frozen[,c("Sample","CTS_frozen")], by = "Sample")

# long format
long_df3 <- df3 %>%
  pivot_longer(
    cols = c(CTS_after, CTS_before, CTS_frozen),
    names_to = "State",
    values_to = "CTS"
  )

# clean labels
long_df3$State <- factor(
  long_df3$State,
  levels = c("CTS_after","CTS_before","CTS_frozen"),
  labels = c("Joint ComBat", "Before ComBat", "Frozen transform")
)

# proportions
library(dplyr)

prop_df3 <- long_df3 %>%
  dplyr::group_by(State, CTS) %>%
  dplyr::summarise(N = n(), .groups = "drop") %>%   # same as count(..., name = "N")
  dplyr::group_by(State) %>%
  dplyr::mutate(Percent = 100 * N / sum(N))


p_hv_prop <- ggplot(prop_df3, aes(State, Percent, fill = CTS)) +
  geom_col(width = 0.6, color = "white", size = 0.6) +
  scale_fill_manual(values = c("1"="royalblue","2"="#B2DF8A","3"="orange")) +
  theme_classic(base_size = 14) +
  labs(title = "Healthy CTS subtype proportions (3 alignment strategies)",
       x = "", y = "Percentage (%)")

ggsave("Figure/CTS_Healthy_Proportion_3methods.svg",
       p_hv_prop, width = 7, height = 5)








## ===============================================================
##  Reorder states & improve aesthetics (Healthy or CAP both apply)
## ===============================================================

# Reorder the factor levels to: Before → Joint → Frozen
long_df3$State <- factor(
  long_df3$State,
  levels = c("Before ComBat", "Joint ComBat", "Frozen transform")
)

# Recalculate proportions after reorder
prop_df3 <- long_df3 %>%
  dplyr::group_by(State, CTS) %>%
  dplyr::summarise(N = n(), .groups = "drop") %>%
  dplyr::group_by(State) %>%
  dplyr::mutate(Percent = 100 * N / sum(N))


## ===============================================================
##  Plot (with adjusted width, spacing, fonts)
## ===============================================================

p_hv_prop <- ggplot(prop_df3, aes(State, Percent, fill = CTS)) +
  
  # wider bars (0.7–0.8 looks good)
  geom_col(width = 0.72, color = "white", size = 0.6) +
  
  scale_fill_manual(values = c("1"="royalblue","2"="#B2DF8A","3"="orange")) +
  
  # consistent font & spacing
  theme_classic(base_size = 14) +
  theme(
    axis.text.x  = element_text(size = 13, angle = 0, hjust = 0.5),
    axis.title.y = element_text(size = 14),
    plot.title   = element_text(size = 17, hjust = 0.5),
    legend.title = element_text(size = 13),
    legend.text  = element_text(size = 12),
    panel.spacing = unit(1.2, "lines")
  ) +
  
  labs(
    title = "Non-infectious CTS subtype proportions (3 alignment strategies)",
    x = "",
    y = "Percentage (%)"
  )

# Save nicely
ggsave("Figure/CTS_Healthy_Proportion_3methods.svg",
       p_hv_prop, width = 6.5, height = 4)
