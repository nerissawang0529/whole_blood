library(dplyr)
library(readr)
library(purrr)
library(broom)

## 0) 先把 value 解析成数值列 value_num（非常关键）
luminex_long2 <- luminex_long %>%
  mutate(
    value_chr = as.character(value),
    # 点作小数、逗号作千分位
    v_dot    = readr::parse_number(value_chr,
                                   locale = readr::locale(decimal_mark = ".", grouping_mark = ",")),
    # 逗号作小数、点作千分位（有些导出会这样）
    v_comma  = readr::parse_number(value_chr,
                                   locale = readr::locale(decimal_mark = ",", grouping_mark = ".")),
    value_num = coalesce(v_dot, v_comma)
  )

# 简单检查
stopifnot(sum(is.finite(luminex_long2$value_num), na.rm = TRUE) > 0)

## 1) 取出 HV 与 CAP 的 LPS–TNF 数值向量
hv_vec <- luminex_long2 %>%
  filter(group == "HV", stimulation == "LPS", marker == "TNF") %>%
  pull(value_num) %>%
  na.omit()

cap_vec <- luminex_long2 %>%
  filter(group == "CAP", stimulation == "LPS", marker == "TNF") %>%
  pull(value_num) %>%
  na.omit()

stopifnot(length(hv_vec) > 0, length(cap_vec) > 0)

hv_mean <- mean(hv_vec)
hv_sd   <- sd(hv_vec)
if (is.na(hv_sd) || hv_sd == 0) stop("HV 的 SD 为 0 或 NA，无法计算 z 分数。")

## 2) 候选 cutoff：使用 CAP 的唯一数值
candidates <- sort(unique(cap_vec))

## A) 均值最接近法：使 mean(CAP<=cutoff) 最接近 hv_mean
res_meanmatch <- tibble(cutoff = candidates) %>%
  mutate(
    n_cap_le    = map_int(cutoff, ~ sum(cap_vec <= .x)),
    mean_cap_le = map_dbl(cutoff, ~ mean(cap_vec[cap_vec <= .x], na.rm = TRUE)),
    diff_mean   = abs(mean_cap_le - hv_mean)
  ) %>%
  filter(n_cap_le >= 10) %>%              # 可按需调整，避免样本过少
  arrange(diff_mean, cutoff) %>%
  slice(1)

cutoff_A <- res_meanmatch$cutoff
z_A      <- (cutoff_A - hv_mean) / hv_sd

print(list(
  method   = "A_mean_match",
  cutoff   = cutoff_A,
  z_score  = z_A,
  hv_mean  = hv_mean,
  hv_sd    = hv_sd,
  n_cap_le = res_meanmatch$n_cap_le,
  mean_cap_le = res_meanmatch$mean_cap_le,
  abs_diff = res_meanmatch$diff_mean
))

## B) 统计学无差异法（Wilcoxon）；选最小使 p>=0.05 的 cutoff
alpha <- 0.05
res_nonsig <- tibble(cutoff = candidates) %>%
  mutate(
    n_cap_le = map_int(cutoff, ~ sum(cap_vec <= .x)),
    p_wilcox = map_dbl(
      cutoff,
      ~ tryCatch(wilcox.test(cap_vec[cap_vec <= .x], hv_vec)$p.value,
                 error = function(e) NA_real_)
    )
  ) %>%
  filter(n_cap_le >= 10, !is.na(p_wilcox)) %>%
  arrange(cutoff) %>%
  filter(p_wilcox >= alpha) %>%
  slice(1)

if (nrow(res_nonsig) > 0) {
  cutoff_B <- res_nonsig$cutoff
  z_B      <- (cutoff_B - hv_mean) / hv_sd
  print(list(
    method   = "B_non_significance",
    cutoff   = cutoff_B,
    z_score  = z_B,
    hv_mean  = hv_mean,
    hv_sd    = hv_sd,
    n_cap_le = res_nonsig$n_cap_le,
    p_wilcox = res_nonsig$p_wilcox
  ))
} else {
  message("方法B：没有找到满足 p>=0.05 的最小 cutoff（可能差异较大或样本量限制）。")
}

## 3) 选择一个 cutoff（这里用方法A的），并给 CAP 病人分组；并回原表
final_cutoff <- cutoff_A

cap_group_df <- luminex_long2 %>%
  filter(group == "CAP", stimulation == "LPS", marker == "TNF") %>%
  transmute(
    ID,
    TNF_LPS_value = value_num,
    TNF_LPS_group = if_else(value_num <= final_cutoff,
                            "Normal_TNF_LPS_CAP", "Elevated_TNF_LPS_CAP")
  )

luminex_long_with_TNFgrp <- luminex_long2 %>%
  left_join(cap_group_df, by = "ID") %>%
  mutate(TNF_LPS_group = if_else(group == "CAP", TNF_LPS_group, NA_character_))

# 简要核对
table(luminex_long_with_TNFgrp$TNF_LPS_group, useNA = "ifany")
final_cutoff
(z_final <- (final_cutoff - hv_mean) / hv_sd)



#box plot
library(dplyr)
library(ggplot2)

cutoff_B <- 8666.78  # 方法B得到的cutoff

# 生成分组
plot_df_B <- luminex_long2 %>%
  filter(stimulation == "LPS", marker == "TNF", group %in% c("HV","CAP")) %>%
  mutate(group_final = case_when(
    group == "HV" ~ "HV",
    group == "CAP" & value_num <= cutoff_B ~ "Normal_CAP",
    group == "CAP" & value_num >  cutoff_B ~ "Elevated_CAP"
  ))

# 绘制 boxplot
p <- ggplot(plot_df_B, aes(x = group_final, y = value_num, fill = group_final)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.8) +
  scale_y_log10() +   # 建议log尺度，便于展示
  labs(title = "LPS-stimulated TNF (Method B cutoff = 8666.78)",
       x = "", y = "TNF (LPS-stimulated)") +
  theme_classic(base_size = 12) +
  scale_fill_manual(values = c("HV" = "#E41A1C",
                               "Normal_CAP" = "#377EB8",
                               "Elevated_CAP" = "#4DAF4A"))

p


plot_df_B$group_final <- factor(plot_df_B$group_final,
                                levels = c("HV", "Normal_CAP", "Elevated_CAP"),
                                labels = c("Non-infection", "Normal LPS TNF", "Elevated LPS TNF"))

p <- ggplot(plot_df_B, aes(x = group_final, y = value_num, fill = group_final)) +
  geom_violin(trim = FALSE, alpha = 0.7, width = 0.8, color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.6, color = "black") +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
  scale_y_log10() +
  labs(title = "LPS-stimulated TNF (cutoff = 8666.78)",
       x = "", y = "TNF (LPS-stimulated)") +
  theme_classic(base_size = 12) +
  scale_fill_manual(values = c("Non-infection" = "#B0B0B0",     # 灰色
                               "Normal LPS TNF" = "#4A90E2",   # 蓝色
                               "Elevated LPS TNF" = "#E67E22")) # 橙色

p
library(dplyr)
library(ggplot2)

# —— 想要的顺序
order_vec <- c("Non-infection", "Normal LPS TNF", "Elevated LPS TNF")

# 1) 先把分组顺序固定
plot_df_B <- plot_df_B %>%
  mutate(group_final = factor(group_final, levels = order_vec))

# 2) 统计样本量，并用“组名(n=xx)”生成标签，保持同样顺序
n_counts <- plot_df_B %>%
  count(group_final) %>%
  mutate(label = paste0(as.character(group_final), "\n(n=", n, ")"))

# 生成一个从组名到标签的映射，按 order_vec 排序
label_map <- n_counts$label[match(order_vec, as.character(n_counts$group_final))]
names(label_map) <- order_vec

# 3) 把标签并回去，并把 label 设为按顺序的因子
plot_df_B <- plot_df_B %>%
  mutate(label = label_map[as.character(group_final)],
         label = factor(label, levels = label_map))   # 关键：这里锁定顺序

# 4) 画图
p <- ggplot(plot_df_B, aes(x = label, y = value_num, fill = group_final)) +
  geom_violin(trim = FALSE, alpha = 0.7, width = 0.8, color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.6, color = "black") +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
  scale_y_log10() +
  labs(title = "Patient stratification based on LPS-stimulated TNF (cutoff = z-score 1.56)",
       x = "", y = "log10 LPS-TNF (pg/ml)") +
  theme_classic(base_size = 14) +
  scale_fill_manual(values = c("Non-infection"   = "#B0B0B0",
                               "Normal LPS TNF" = "#4A90E2",
                               "Elevated LPS TNF" = "#E67E22")) +
  guides(fill = guide_legend(title = "group_final"))

p
