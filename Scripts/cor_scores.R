rm(list = ls ())

## =========================
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
               tableone, kmed, survival, survminer, skimr, lattice, reshape2, grDevices, ggrepel)


SRSq <- import("Original_data/studydata_with_SRSq7.csv") 
SRSq_score <- SRSq[,c("EB_id","SRSq_7")]
CTS <- import("Original_data/studydata_with_cts.csv")   
CTS_score <- CTS[,c("EB_id","CTS")]
MS1 <- import("Original_data/studydata_with_ms1.csv") 
MS1_score <- MS1[,c("EB_id","MS1")]
il1r2 <- import("Original_data/studydata_with_il1r2.csv")
il1r2_score <- il1r2[,c("EB_id","IL1R2")]

merged_scores <- list(SRSq_score, CTS_score, MS1_score, il1r2_score) |>
  purrr::reduce(dplyr::full_join, by = "EB_id")

# Packages
library(dplyr)
library(ggplot2)
library(patchwork)   # 用于拼图

# === 数据准备 ===
df <- merged_scores %>%
  transmute(
    EB_id,
    SRSq = SRSq_7,
    MS1  = MS1,
    CTS  = factor(CTS, levels = sort(unique(CTS)), ordered = TRUE),
    CTS_num = as.numeric(CTS)   # 把有序因子转为 1/2/3 用于 Spearman
  )

# === 计算三组相关 + FDR ===
cor_list <- list(
  c("MS1","SRSq"),
  c("MS1","CTS_num"),
  c("SRSq","CTS_num")
)

cors <- lapply(cor_list, function(v){
  x <- v[1]; y <- v[2]
  xx <- df[[x]]; yy <- df[[y]]
  keep <- complete.cases(xx, yy)
  tst <- suppressWarnings(cor.test(xx[keep], yy[keep], method = "spearman"))
  data.frame(x = x, y = y, rho = unname(tst$estimate), p = tst$p.value)
}) %>% bind_rows()

cors$FDR <- p.adjust(cors$p, method = "BH")

# 取数值方便标注
rho_ms1_srsq  <- cors$rho[cors$x=="MS1"  & cors$y=="SRSq"]
fdr_ms1_srsq  <- cors$FDR[cors$x=="MS1"  & cors$y=="SRSq"]
rho_ms1_cts   <- cors$rho[cors$x=="MS1"  & cors$y=="CTS_num"]
fdr_ms1_cts   <- cors$FDR[cors$x=="MS1"  & cors$y=="CTS_num"]
rho_srsq_cts  <- cors$rho[cors$x=="SRSq" & cors$y=="CTS_num"]
fdr_srsq_cts  <- cors$FDR[cors$x=="SRSq" & cors$y=="CTS_num"]

# === 作图 ===

# 1) MS1 vs SRSq（散点 + 拟合线）
p1 <- ggplot(df, aes(x = SRSq, y = MS1)) +
  geom_point(alpha = .6) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = "SRSq", y = "MS1 (%)",
       title = sprintf("MS1 vs SRSq  (Spearman ρ=%.2f; FDR=%.3g)",
                       rho_ms1_srsq, fdr_ms1_srsq)) +
  theme_classic()

# 2) MS1 vs CTS（箱线/小提琴 + ρ/FDR 以 CTS_num 计）
p2 <- ggplot(df, aes(x = CTS, y = MS1)) +
  geom_violin(trim = FALSE, alpha = .7) +
  geom_boxplot(width = .15, outlier.shape = NA) +
  labs(x = "CTS (ordered)", y = "MS1 (%)",
       title = sprintf("MS1 across CTS  (Spearman ρ=%.2f; FDR=%.3g)",
                       rho_ms1_cts, fdr_ms1_cts)) +
  theme_classic()

# 3) SRSq vs CTS（箱线/小提琴 + ρ/FDR）
p3 <- ggplot(df, aes(x = CTS, y = SRSq)) +
  geom_violin(trim = FALSE, alpha = .7) +
  geom_boxplot(width = .15, outlier.shape = NA) +
  labs(x = "CTS (ordered)", y = "SRSq",
       title = sprintf("SRSq across CTS  (Spearman ρ=%.2f; FDR=%.3g)",
                       rho_srsq_cts, fdr_srsq_cts)) +
  theme_classic()

# === 拼成“三角相关”面板 ===
p4 <- (p1 | p2 | p3)

ggsave("Figure/cor_score.svg", p4, width = 15, height = 6, dpi = 300)






# correlation MS1 vs IL1R2 
cor_s <- cor.test(merged_scores$MS1, merged_scores$IL1R2, method = "spearman")

p6 <- ggplot(merged_scores, aes(x = MS1, y = IL1R2)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", color = "#0072B2") +
  labs(
    title = sprintf("MS1 vs IL1R2 (Spearman ρ = %.2f, p = %.2g)",
                    cor_s$estimate, cor_s$p.value),
    x = "MS1 score", y = "IL1R2 (%)"
  ) +
  theme_classic(base_size = 12)

ggsave("Figure/cor_score.svg", p6, width = 15, height = 6, dpi = 300)








# Make sure numeric
merged_scores <- merged_scores %>%
  mutate(
    IL1R2 = as.numeric(IL1R2),
    MS1   = as.numeric(MS1)
  )

# Keep only IL1R2 > 0 for correlation/plot
df_pos <- merged_scores %>%
  filter(!is.na(IL1R2), IL1R2 > 0, !is.na(MS1))

# Spearman among non-zeros
cor_s <- cor.test(df_pos$MS1, df_pos$IL1R2, method = "spearman", exact = FALSE)

# Plot
p6_nonzero <- ggplot(df_pos, aes(x = MS1, y = IL1R2)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", se = TRUE) +
  labs(
    title = sprintf("MS1 vs IL1R2 (non-zero only) — Spearman \u03C1 = %.2f, p = %.2g",
                    cor_s$estimate, cor_s$p.value),
    x = "MS1 score", y = "IL1R2 (%)"
  ) +
  theme_classic(base_size = 12)

ggsave("Figure/cor_score_nonzero.svg", p6_nonzero, width = 8, height = 6, dpi = 300)
