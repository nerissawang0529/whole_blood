rm(list = ls())
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
               tableone, kmed, survival, survminer, skimr, lattice, reshape2, grDevices, ggrepel)
install.packages("conflicted")  # 只需安装一次
library(conflicted)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::mutate)
conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::arrange)
conflicts_prefer(dplyr::lag)     # 需要的话


#ready data
luminex_long <- read.csv("Original_data/luminex_long_1.csv")

DIP_continous <- read.csv("Original_data/DIP_predictions_20250606_continous.csv")

DIP_stage <- read.csv("Original_data/DIP_predictions_20250606stage.csv")

## ============== Packages ==============
suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

## ============== 0) 规范列名（稳妥） ==============
# 目标：
# DIP_continous:  ID, cDIP
# DIP_stage:      ID, DIP, DIP1_Prob, DIP2_Prob, DIP3_Prob

# -- DIP_continous --
nm1 <- names(DIP_continous)
col_ID1  <- nm1[grepl("^ID$", nm1, ignore.case = TRUE)][1]
col_cDIP <- nm1[grepl("^(c?DIP)$|^cDIP$", nm1, ignore.case = TRUE)][1]
DIP_continous_std <- DIP_continous %>%
  dplyr::rename(ID = !!col_ID1, cDIP = !!col_cDIP)

# -- DIP_stage --
nm2 <- names(DIP_stage)
col_ID2 <- nm2[grepl("^ID$", nm2, ignore.case = TRUE)][1]
col_DIP <- nm2[grepl("^DIP$", nm2, ignore.case = TRUE)][1]
col_P1  <- nm2[grepl("^DIP1[_\\.]?Prob$", nm2, ignore.case = TRUE)][1]
col_P2  <- nm2[grepl("^DIP2[_\\.]?Prob$", nm2, ignore.case = TRUE)][1]
col_P3  <- nm2[grepl("^DIP3[_\\.]?Prob$", nm2, ignore.case = TRUE)][1]

DIP_stage_std <- DIP_stage %>%
  dplyr::rename(
    ID        = !!col_ID2,
    DIP       = !!col_DIP,
    DIP1_Prob = !!col_P1,
    DIP2_Prob = !!col_P2,
    DIP3_Prob = !!col_P3
  )

## ============== 1) 生成最终 DIP 分组 + 合并（仅保留DIP里出现过的ID） ==============
DIP_stage2 <- DIP_stage_std %>%
  mutate(
    DIP_confidence  = pmax(DIP1_Prob, DIP2_Prob, DIP3_Prob, na.rm = TRUE),
    DIP_top_by_prob = c("DIP1","DIP2","DIP3")[
      max.col(cbind(DIP1_Prob, DIP2_Prob, DIP3_Prob), ties.method = "first")
    ],
    DIP_stage_final = ifelse(is.na(DIP) | !DIP %in% c("DIP1","DIP2","DIP3"),
                             DIP_top_by_prob, DIP)
  ) %>%
  dplyr::select(ID, DIP_stage_final, DIP_confidence,
                DIP1_Prob, DIP2_Prob, DIP3_Prob)

# 只保留在“连续分数”和“分组”两张表都出现的 ID（避免 NA）
dip_info <- DIP_continous_std %>%
  inner_join(DIP_stage2, by = "ID")

## ============== 2) 取 CAP + LPS + TNF，并按 TNF 升序；只保留 DIP 有ID ==============
tnf_cap_lps <- luminex_long %>%
  filter(group == "CAP",
         stimulation == "LPS") %>%
  # 若你的 TNF 是精确名字（如 "TNF_1"），把下一行改为：filter(marker == "TNF_1")
  filter(str_detect(marker, regex("TNF", ignore_case = TRUE))) %>%
  mutate(
    value_num = readr::parse_number(as.character(value),
                                    locale = readr::locale(decimal_mark = ".", grouping_mark = ","))
  ) %>%
  # 关键：只保留在 DIP 里出现过的 ID（inner_join）
  inner_join(dip_info, by = "ID") %>%
  arrange(value_num) %>%
  mutate(
    rank = row_number(),
    TNF_log10 = log10(value_num),
    DIP_stage_final = factor(DIP_stage_final, levels = c("DIP1","DIP2","DIP3"))
  ) %>%
  tidyr::drop_na(DIP_stage_final)   # 仍防守性去掉极少数缺失

## ============== 3) 两条热图（上：cDIP；下：DIP stage） ==============
# 视觉改良：统一极简风、左对齐标题、合并图例、去掉厚重边框
# 连续色带用更平滑的渐变（不刺眼），离散色用柔和三色

p_heat_cdip <- ggplot(tnf_cap_lps, aes(x = rank, y = "cDIP (continuous)", fill = cDIP)) +
  geom_tile(height = 0.98) +
  scale_fill_gradientn(
    colours = c("#f7fbff","#c6dbef","#6baed6","#2171b5","#08306b"),
    name = "cDIP"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = NULL, title = "LPS-TNF-ordered heatmaps (CAP patients)") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0)
  )

p_heat_stage <- ggplot(tnf_cap_lps, aes(x = rank, y = "DIP (stage)", fill = DIP_stage_final)) +
  geom_tile(height = 0.98) +
  scale_fill_manual(values = c(DIP1 = "#4C78A8", DIP2 = "#F58518", DIP3 = "#54A24B"),
                    name = "DIP stage") +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  labs(x = "Patients sorted by TNF (ascending)", y = NULL) +
  theme_classic(base_size = 12) +
  theme(legend.position = "right")

p_heat <- p_heat_cdip / p_heat_stage + plot_layout(heights = c(1, 1), guides = "collect")
print(p_heat)

## ============== 4) 相关性图（美化） ==============
cor_s <- cor.test(tnf_cap_lps$TNF_log10, tnf_cap_lps$cDIP, method = "spearman")

ann_txt <- sprintf("Spearman \u03c1 = %.2f (p = %.2g)",
                   unname(cor_s$estimate), cor_s$p.value)

p_scatter <- ggplot(tnf_cap_lps, aes(x = TNF_log10, y = cDIP)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.7,
              color = "black", fill = "grey80") +
  annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.5,
           label = ann_txt, size = 3.6) +
  labs(x = "TNF (log10) • LPS • CAP", y = "cDIP (continuous)",
       title = "Association of TNF and cDIP (Spearman)") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0))

print(p_scatter)
