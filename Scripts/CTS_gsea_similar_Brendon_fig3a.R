## 
rm(list = ls())
pkgs <- c("tidyverse","limma","AnnotationDbi","org.Hs.eg.db","msigdbr","fgsea","ggtext")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
invisible(lapply(pkgs, library, character.only = TRUE))


# --- Load data ---
CTS_predictions <- read.csv("Original_data/CTS_predictions.csv", check.names = FALSE)
expr_logcpm_df  <- read.csv("Original_data/expr_logcpm.csv", row.names = 1, check.names = FALSE)
expr <- as.matrix(expr_logcpm_df)
expr <- expr[, CTS_predictions$Sample, drop = FALSE]

cts <- factor(CTS_predictions$CTS, levels = c(1,2,3))
design <- model.matrix(~0 + cts)
colnames(design) <- c("CTS1","CTS2","CTS3")


## ===== 1. 对齐样本 =====
keep <- intersect(CTS_predictions$Sample, colnames(expr_logcpm_df))
cts_df <- CTS_predictions %>%
  filter(Sample %in% keep) %>%
  mutate(CTS = factor(CTS, levels = c(1,2,3),
                      labels = c("CTS1","CTS2","CTS3")))
expr <- as.matrix(expr_logcpm_df[, keep, drop = FALSE])
rownames(expr) <- sub("\\..*$","", rownames(expr))   # 去掉版本号
storage.mode(expr) <- "double"

## ===== 2. Ensembl → Entrez 映射 =====
entrez_map <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                    keys = unique(rownames(expr)),
                                    keytype = "ENSEMBL",
                                    column  = "ENTREZID",
                                    multiVals = "first")

## ===== 3. limma 差异分析（CTS1/2/3 vs others）=====
design <- model.matrix(~0 + CTS, data = cts_df)
colnames(design) <- levels(cts_df$CTS)

contr <- makeContrasts(
  CTS1_vs_others = CTS1 - (CTS2 + CTS3)/2,
  CTS2_vs_others = CTS2 - (CTS1 + CTS3)/2,
  CTS3_vs_others = CTS3 - (CTS1 + CTS2)/2,
  levels = design)

fit <- lmFit(expr, design)
fit2 <- contrasts.fit(fit, contr)
fit3 <- eBayes(fit2)
tmat <- fit3$t     # 每个对比的t统计量矩阵

## ===== 4. Hallmark基因集 =====
hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, entrez_gene)
hallmark_sets <- split(hallmark$entrez_gene, hallmark$gs_name)

## ===== 5. 构建geneList并运行GSEA =====
make_geneList <- function(t_vec, ens_ids, entrez_map) {
  tibble(ENSEMBL = ens_ids, t = as.numeric(t_vec)) %>%
    mutate(ENTREZ = entrez_map[ENSEMBL]) %>%
    filter(!is.na(ENTREZ)) %>%
    group_by(ENTREZ) %>%
    slice_max(order_by = abs(t), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    { setNames(.$t, .$ENTREZ) } %>%
    sort(decreasing = TRUE)
}

set.seed(1)
gsea_all <- lapply(colnames(tmat), function(cn) {
  gl <- make_geneList(tmat[, cn], rownames(tmat), entrez_map)   # 已经是命名的数值向量、降序
  fg <- fgsea(pathways = hallmark_sets,
              stats    = gl,
              minSize  = 10,
              maxSize  = 500,
              scoreType = "std") %>%                             # 可省略，保持默认
    arrange(padj) %>% 
    mutate(contrast = cn)
  fg
}) %>% bind_rows()


## ===== 6. 整理数据用于绘图 =====
# 把 fgsea 输出列 'pathway' 改名为 'gs_name'
gsea_all <- gsea_all %>% dplyr::rename(gs_name = pathway)

plot_df <- gsea_all %>%
  mutate(FDRneglog10 = -log10(padj),
         sig = padj < 0.05,
         gs_name = str_remove(gs_name, "^HALLMARK_") %>% str_replace_all("_", " ")) %>%
  group_by(contrast) %>%
  arrange(desc(sig), desc(FDRneglog10)) %>%
  slice_head(n = 20) %>%
  ungroup() %>%
  mutate(gs_name = factor(gs_name, levels = rev(unique(gs_name))),
         contrast = recode(contrast,
                           CTS1_vs_others = "CTS1",
                           CTS2_vs_others = "CTS2",
                           CTS3_vs_others = "CTS3"))

## ===== 7. 绘制Fig 3a风格图 =====
ggplot(plot_df, aes(x = FDRneglog10, y = gs_name)) +
  geom_col(fill = "grey50", width = 0.7) +
  geom_point(aes(size = sig), color = "black") +
  scale_size_manual(values = c(`TRUE` = 2.8, `FALSE` = 1.5), guide = "none") +
  facet_wrap(~ contrast, ncol = 3, scales = "free_y") +
  labs(x = expression(-log[10]("FDR q")), y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_markdown(size = 8),
    axis.text.x = element_text(size = 10),
    strip.text = element_text(size = 12, face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

## ===== 8. 保存结果 =====
ggsave("Figure3a_CAP_CTS_GSEA.svg", width = 10, height = 6, dpi = 300)

library(tidyverse)
library(ggtext)

# ---- 0) 清洗 + 只留显著 ----
df <- gsea_all %>%
  mutate(
    FDRneglog10 = -log10(padj),
    pathway     = stringr::str_remove(gs_name, "^HALLMARK_") %>% stringr::str_replace_all("_"," "),
    contrast    = dplyr::recode(contrast,
                                CTS1_vs_others = "CTS1",
                                CTS2_vs_others = "CTS2",
                                CTS3_vs_others = "CTS3")
  ) %>%
  filter(padj < 0.05)

# ---- 1) 每个 CTS 取前 n 个通路（按显著性或 |NES| 任选其一）----
n_keep <- 5
df_top <- df %>%
  dplyr::group_by(contrast) %>%
  dplyr::arrange(desc(FDRneglog10), .by_group = TRUE) %>%   # 若想按 |NES|：arrange(desc(abs(NES)), .by_group=TRUE)
  slice_head(n = n_keep) %>%
  ungroup()

# ---- 2) 生成唯一的 y 轴键，避免同名通路在不同 CTS 上“叠加” ----
# 显示仍用纯 pathway 名称；排序 = 先 CTS1 块、再 CTS2、再 CTS3；块内按显著性从小到大（画出来自下而上递增）
order_levels <- df_top %>%
  mutate(contrast = factor(contrast, levels = c("CTS1","CTS2","CTS3"))) %>%
  arrange(contrast, FDRneglog10) %>%
  transmute(key = sprintf("%s__%s", pathway, contrast)) %>%
  pull(key)

df_top <- df_top %>%
  mutate(
    key      = sprintf("%s__%s", pathway, contrast),           # 唯一键，避免叠加
    key      = factor(key, levels = unique(order_levels)),
    label_y  = pathway                                         # 仅用于显示的标签
  )

# ---- 3) 颜色（与论文一致：蓝=CTS1，绿=CTS2，橙=CTS3）----
cts_cols <- c("CTS1"="#355E9A","CTS2"="#9ACD66","CTS3"="#F5A623")

# ---- 4) 作图：单面板、每条仅一个 CTS、颜色=CTS、横轴= -log10(FDR q) ----
p <- ggplot(df_top, aes(x = FDRneglog10, y = key, fill = contrast)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = cts_cols, name = NULL) +
  scale_y_discrete(labels = setNames(df_top$label_y, df_top$key)) +
  labs(x = expression("FDR q ("*-log[10]*")"), y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    axis.text.y  = ggtext::element_markdown(size = 10),
    axis.text.x  = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank()
  )

# 如果你想让横轴和论文一样大约到 5，可以加上限制（可选）：
# p <- p + scale_x_continuous(limits = c(0, 5))

print(p)
ggsave("Figure3a_singlepanel_CAP_GSEA_final.svg", p, width = 6.5, height = 6.5, dpi = 300)


## 放在 df_top 创建完成之后（作图之前）

# 1) 明确 CTS 顺序
df_top <- df_top %>%
  mutate(contrast = factor(contrast, levels = c("CTS1","CTS2","CTS3")))

# 2) 生成按 CTS1→CTS2→CTS3、块内按显著性排序的 y 轴顺序
order_levels <- df_top %>%
  arrange(contrast, desc(FDRneglog10)) %>%                 # 块顺序 + 块内从显著到不显著
  transmute(key = sprintf("%s__%s", pathway, contrast)) %>%
  pull(key)

# 3) 用这个顺序作为 y 轴因子水平（CTS1 在最上方：用 rev；若想 CTS1 在最下方：去掉 rev）
df_top <- df_top %>%
  mutate(key = factor(sprintf("%s__%s", pathway, contrast),
                      levels = rev(order_levels)))         # ← CTS1 顶部

ggplot(df_top, aes(x = FDRneglog10, y = key, fill = contrast)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.2) +
  scale_y_discrete(labels = setNames(df_top$pathway, df_top$key)) +
  scale_fill_manual(values = c("CTS1"="#355E9A","CTS2"="#9ACD66","CTS3"="#F5A623"), name = NULL) +
  labs(x = expression("FDR q ("*-log[10]*")"), y = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        axis.text.y = ggtext::element_markdown(size = 10),
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank())


library(stringi)
library(stringr)
## ---- 生成用于显示的标签（句首大写）----
df_top <- df_top %>%
  mutate(
    # 先全小写再句首大写 -> “Protein secretion” 风格
    label_y = str_to_sentence(str_to_lower(pathway))
  )

## ---- 恢复常见缩写/专有名大小写（可按需增减）----
fix_map <- c(
  "Il6" = "IL-6",
  "Jak" = "JAK",
  "Stat" = "STAT",
  "Stat3" = "STAT3",
  "Dna" = "DNA",
  "Rna" = "RNA",
  "Mtorc1" = "mTORC1",
  "Mtorc2" = "mTORC2",
  "Kras" = "KRAS",
  "Myc"  = "MYC",
  "Tnfa" = "TNFα",
  "Wnt/β-catenin" = "Wnt/β-catenin",
  "Wnt beta catenin" = "Wnt β-catenin",
  "Nfkb" = "NF-κB"
)

for (k in names(fix_map)) {
  df_top$label_y <- str_replace_all(df_top$label_y, fixed(k), fix_map[[k]])
}

## ---- 用新标签画图 ----
ggplot(df_top, aes(x = FDRneglog10, y = key, fill = contrast)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.2) +
  scale_y_discrete(labels = setNames(df_top$label_y, df_top$key)) +
  scale_fill_manual(values = c("CTS1"="#355E9A","CTS2"="#9ACD66","CTS3"="#F5A623"), name = NULL) +
  labs(x = expression("FDR q ("*-log[10]*")"), y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    axis.text.y  = ggtext::element_markdown(size = 10),
    axis.text.x  = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank()
  )


cts_cols <- c("CTS1"="#355E9A","CTS2"="#9ACD66","CTS3"="#F5A623")


df_best <- gsea_all %>%
  as.data.frame() %>%
  mutate(
    pathway  = stringr::str_remove(gs_name, "^HALLMARK_") %>% stringr::str_replace_all("_"," "),
    contrast = dplyr::recode(contrast,
                             CTS1_vs_others="CTS1",
                             CTS2_vs_others="CTS2",
                             CTS3_vs_others="CTS3"),
    FDRneglog10 = -log10(padj)
  ) %>%
  # 只保留有有效q值且显著的
  dplyr::filter(!is.na(padj), padj < 0.05) %>%
  # 每条通路挑选 q 最小的那一行；若并列，再用 |NES| 大的优先作为决胜
  dplyr::group_by(pathway) %>%
  dplyr::arrange(padj, dplyr::desc(abs(NES)), .by_group = TRUE) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()


ggplot(df_best, aes(x = FDRneglog10, y = pathway, fill = contrast)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.25) +
  scale_fill_manual(values = cts_cols, name = NULL) +
  labs(x = expression("FDR q ("*-log[10]*")"), y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    axis.text.y  = element_text(size = 10),
    axis.text.x  = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank()
  )


# —— 排序与标签整理（保持你已有的 df_best 计算不变）——
df_best2 <- df_best %>%
  dplyr::mutate(
    contrast = factor(contrast, levels = c("CTS1","CTS2","CTS3")),  # 块顺序：CTS1→CTS2→CTS3
    pathway_label = stringr::str_to_sentence(pathway)               # 仅首字母大写，其余小写
  ) %>%
  dplyr::arrange(contrast, dplyr::desc(FDRneglog10)) %>%            # 先按 CTS 分块，再块内按显著性降序
  dplyr::mutate(
    # 让“当前行顺序”显示在 y 轴从上到下（factor 的 levels 从后到前显示）
    pathway_ord = factor(pathway_label, levels = rev(pathway_label))
  )

# —— 画图 ——（cts_cols 用你之前定义的配色：CTS1 蓝、CTS2 绿、CTS3 橙）
ggplot(df_best2, aes(x = FDRneglog10, y = pathway_ord, fill = contrast)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.25) +
  scale_fill_manual(values = cts_cols, name = NULL) +
  labs(x = expression("FDR q ("*-log[10]*")"), y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    axis.text.y  = element_text(size = 10),
    axis.text.x  = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank()
  )
## ===========================================
##  EXPORT FINAL GSEA FIGURE  (df_best2)
## ===========================================
dir.create("Figure", showWarnings = FALSE)

p_final <- ggplot(df_best2, aes(x = FDRneglog10, y = pathway_ord, fill = contrast)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.25) +
  scale_fill_manual(values = c("CTS1"="#355E9A",
                               "CTS2"="#9ACD66",
                               "CTS3"="#F5A623"), name = NULL) +
  labs(x = expression("FDR q ("*-log[10]*")"), y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    axis.text.y  = element_text(size = 10),
    axis.text.x  = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank()
  )

# Preview
print(p_final)

# ==== Save as SVG ====
ggsave("Figure/CTS_CAP_GSEA_final.svg",
       p_final, width = 6.5, height = 7, dpi = 300, device = "svg")

# ==== (Optional) Save as PDF ====
ggsave("Figure/CTS_CAP_GSEA_final.pdf",
       p_final, width = 6.5, height = 7)

