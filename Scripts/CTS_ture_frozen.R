#need to run the code 'CTS_all_RNA_makescore_table1' and 'CTS_CAP_before_after_combat'first

## ===============================================================
## 1b. Build a simple FROZEN batch transform from sepsis only
##     (location + scale per gene, ComBat-style but very simple)
## ===============================================================

# Use the same common genes as in your RF code
common_genes <- intersect(rownames(exp_core_g), rownames(expr_logcpm))

# Sepsis reference (already ComBat-normalized by Brendon)
core_mat <- exp_core_g[common_genes, , drop = FALSE]

# Gene-wise mean and SD in sepsis (this defines the FROZEN "sepsis space")
mu_core <- rowMeans(core_mat)
sd_core <- apply(core_mat, 1, sd)
sd_core[sd_core == 0 | is.na(sd_core)] <- 1  # avoid zero/NA


# CAP raw logCPM, restricted to the same genes
cap_mat_raw <- expr_logcpm[common_genes, , drop = FALSE]

# Gene-wise mean and SD in CAP (current batch)
mu_cap <- rowMeans(cap_mat_raw)
sd_cap <- apply(cap_mat_raw, 1, sd)
sd_cap[sd_cap == 0 | is.na(sd_cap)] <- 1

# ---- FROZEN TRANSFORM: put CAP onto the sepsis (core) scale ----
# 1) center CAP by its own mean
cap_frozen <- sweep(cap_mat_raw, 1, mu_cap, "-")
# 2) standardize by its own SD
cap_frozen <- sweep(cap_frozen, 1, sd_cap, "/")
# 3) rescale to sepsis SD
cap_frozen <- sweep(cap_frozen, 1, sd_core, "*")
# 4) shift to sepsis mean
cap_frozen <- sweep(cap_frozen, 1, mu_core, "+")


## ===============================================================
## 1c. Train RF on sepsis (same as before) and predict on FROZEN CAP
## ===============================================================

y_core <- factor(core_samples[colnames(core_mat), "CTS"])

set.seed(123)
rf_model_frozen <- randomForest(x = t(core_mat), y = y_core)

pred_frozen <- data.frame(
  Sample   = colnames(cap_frozen),
  CTS_frozen = predict(rf_model_frozen, newdata = t(cap_frozen))
)

pred_noCombat <- data.frame(Sample = colnames(cap_mat),
                            CTS_noCombat = predict(rf_model, newdata = t(cap_mat)))

pred_combat <- res$predictions[, c("Sample","CTS")]
names(pred_combat) <- c("Sample","CTS_ComBat")

cmp <- inner_join(pred_noCombat, pred_combat, by = "Sample")
cmp <- cmp %>%
  inner_join(pred_frozen, by = "Sample")
cmp_long3 <- cmp %>%
  pivot_longer(cols = c(CTS_noCombat, CTS_ComBat, CTS_frozen),
               names_to = "Version",
               values_to = "CTS")

cmp_long3$Version <- factor(
  cmp_long3$Version,
  levels = c("CTS_noCombat", "CTS_ComBat", "CTS_frozen"),
  labels = c("Before ComBat", "Joint ComBat", "Frozen transform")
)

p_prop3 <- ggplot(cmp_long3, aes(Version, fill = CTS)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("1"="royalblue","2"="#B2DF8A","3"="orange")) +
  theme_classic(base_size = 12) +
  labs(title = "CAP CTS subtype proportions (alignment strategies)",
       y = "Proportion", x = NULL)

ggsave("Results/CTS_Proportion_Before_After_Frozen_ComBat.svg", p_prop3, width = 6.5, height = 4)
ggsave("Results/CTS_Proportion_Before_After_Frozen_ComBat.pdf", p_prop3, width = 6.5, height = 4)



## ===== 计算每个 Version × CTS 的比例 =====
library(dplyr)

prop_cmp3 <- cmp_long3 %>%
  dplyr::group_by(Version, CTS) %>%
  dplyr::summarise(N = n(), .groups = "drop") %>%
  dplyr::group_by(Version) %>%
  dplyr::mutate(Percent = 100 * N / sum(N))

## ===== 画带百分比标注的堆叠条形图 =====
p_prop3 <- ggplot(prop_cmp3, aes(Version, Percent, fill = CTS)) +
  geom_col(width = 0.9, color = "white", size = 0.6) +          # 和 HV 那张一样的宽度风格
  geom_text(
    aes(label = sprintf("%.0f%%", Percent)),                    # 显示整数百分比，可改成"%.1f%%"
    position = position_stack(vjust = 0.5),                     # 放在每个色块中间
    color = "white",
    size  = 4
  ) +
  scale_y_continuous(
    name   = "Proportion",                                      # 或者改成 "Percentage (%)"
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_fill_manual(values = c("1"="royalblue","2"="#B2DF8A","3"="orange")) +
  theme_classic(base_size = 12) +
  labs(
    title = "CAP CTS subtype proportions (alignment strategies)",
    x = NULL
  )

ggsave("Results/CTS_Proportion_Before_After_Frozen_ComBat.svg",
       p_prop3, width = 6.5, height = 4)
ggsave("Results/CTS_Proportion_Before_After_Frozen_ComBat.pdf",
       p_prop3, width = 6.5, height = 4)
