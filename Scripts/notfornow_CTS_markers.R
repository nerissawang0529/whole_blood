





# make the primery figure#####
library(ggplot2)
library(dplyr)

# 1) Set order + short labels
wb2 <- wb %>%
  mutate(
    stimulation = factor(
      stimulation,
      levels = c("M","KP","PP","LPS"),
      labels = c("Medium","Kleb.","Pneumo.","LPS")
    )
  )

# 2) Plot with dodge + jitter-dodge and filled boxes (clearer than outlines)
p_primary <- ggplot(wb2, aes(x = stimulation, y = value, fill = CTS)) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.6, alpha = 0.55, colour = "#4a4a4a",
    outlier.shape = NA
  ) +
  geom_point(
    aes(color = CTS),
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
    size = 1.4, alpha = 0.65
  ) +
  scale_y_continuous(trans = "log10") +
  facet_wrap(~ marker, scales = "free_y", ncol = 3) +
  labs(x = NULL, y = "Cytokine (log10 scale)",
       title = "Whole-blood stimulation: primary readouts by CTS") +
  theme_classic(base_size = 12) +
  theme(strip.text = element_text(face = "bold"),
        legend.position = "top") +
  guides(color = "none")   # keep one legend (from fill)

p_primary
library(dplyr)
library(rstatix)
library(ggpubr)

wb3 <- wb %>%
  mutate(
    stimulation = factor(stimulation,
                         levels = c("M","KP","PP","LPS"),
                         labels = c("Medium","Kleb.","Pneumo.","LPS")
    ),
    CTS = factor(CTS)   # make sure it's a factor: 1,2,3
  )

my_comps <- list(c("1","2"), c("1","3"), c("2","3"))

stat_df <- wb3 %>%
  group_by(marker, stimulation) %>%
  group_modify(~{
    d <- .x
    y_max <- max(d$value, na.rm = TRUE)
    
    wilcox_test(d, value ~ CTS, comparisons = my_comps, p.adjust.method = "BH") %>%
      add_significance(p.col = "p.adj") %>%
      mutate(
        pair = paste(group1, group2, sep = "-"),
        rank = match(pair, c("1-2","1-3","2-3")),
        y.position = y_max * (1.15 + 0.08 * (rank - 1))
      )
  }) %>%
  ungroup()

p_sig <- ggplot(wb3, aes(x = CTS, y = value, fill = CTS)) +
  geom_boxplot(width = 0.65, alpha = 0.55, colour = "#4a4a4a", outlier.shape = NA) +
  geom_jitter(aes(color = CTS), width = 0.12, size = 1.2, alpha = 0.6) +
  scale_y_continuous(trans = "log10") +
  facet_grid(marker ~ stimulation, scales = "free_y") +
  labs(x = NULL, y = "Cytokine (log10 scale)",
       title = "Whole-blood stimulation by CTS (with BH-adjusted significance)") +
  theme_classic(base_size = 12) +
  theme(strip.text = element_text(face = "bold"),
        legend.position = "top") +
  guides(color = "none") +
  ggpubr::stat_pvalue_manual(
    stat_df,
    label = "p.adj.signif",            # shows ****/***/**/*
    hide.ns = TRUE,                    # drop 'ns' to reduce clutter
    tip.length = 0.01,
    bracket.size = 0.3,
    size = 3
  )

p_sig#####


## ===== Pairwise contrasts with robust direction (first − second) =====
library(dplyr); library(forcats); library(rstatix); library(purrr); library(ggplot2)

cts_lvls <- sort(unique(as.character(wb_clean[[group_var]])))

# 我们定义列标题顺序（第一组在前）
contrast_pairs <- list(
  c("1","2"),  # CTS1 vs CTS2
  c("1","3"),  # CTS1 vs CTS3
  c("2","3")   # CTS2 vs CTS3
)
contrast_df <- as.data.frame(do.call(rbind, contrast_pairs))
names(contrast_df) <- c("first","second")

effsize_df <- purrr::pmap_dfr(contrast_df, function(first, second) {
  
  dat <- wb_clean %>%
    filter(.data[[group_var]] %in% c(first, second)) %>%
    mutate(!!group_var := fct_relevel(.data[[group_var]], first, second))
  
  # 1) 先按 rstatix 计算（它默认第二组−第一组的方向）
  es <- dat %>%
    group_by(stimulation, marker) %>%
    cohens_d(reformulate(group_var, response = "value_num"),
             var.equal = FALSE, hedges.correction = TRUE, ci = TRUE) %>%
    ungroup() %>%
    rename(g_raw = effsize, lo_raw = conf.low, hi_raw = conf.high)
  
  # 2) 用“第一组−第二组”的均值差，来校准方向
  md <- dat %>%
    group_by(stimulation, marker, !!sym(group_var)) %>%
    summarise(mu = mean(value_num, na.rm = TRUE), .groups = "drop_last") %>%
    summarise(delta = mu[.data[[group_var]] == first] - mu[.data[[group_var]] == second],
              .groups = "drop")
  
  out <- es %>%
    left_join(md, by = c("stimulation","marker")) %>%
    mutate(
      flip = ifelse(sign(delta) == 0 | is.na(delta), 1,
                    ifelse(sign(delta) == sign(g_raw), 1, -1)),
      g        = g_raw  * flip,          # 方向校正
      conf.low = lo_raw * flip,
      conf.high= hi_raw * flip,
      contrast_lab = sprintf("CTS%s vs CTS%s", first, second)
    ) %>%
    select(stimulation, marker, g, conf.low, conf.high, contrast_lab)
  
  # 3) p 值与显著性
  tt <- dat %>%
    group_by(stimulation, marker) %>%
    t_test(reformulate(group_var, response = "value_num"), var.equal = FALSE) %>%
    ungroup() %>%
    mutate(p_BH = p.adjust(p, method = "BH")) %>%
    select(stimulation, marker, p_BH)
  
  out %>%
    left_join(tt, by = c("stimulation","marker")) %>%
    mutate(
      p_BH.signif = case_when(
        is.na(p_BH) ~ "ns",
        p_BH < 1e-4 ~ "****",
        p_BH < 1e-3 ~ "***",
        p_BH < 1e-2 ~ "**",
        p_BH < 5e-2 ~ "*",
        TRUE ~ "ns"
      ),
      dir = case_when(g > 0 ~ "UP", g < 0 ~ "DOWN", TRUE ~ NA_character_),
      signif_dir = ifelse(p_BH.signif == "ns" | is.na(dir), "ns", paste(p_BH.signif, dir))
    )
})

## ===== 标签映射 =====
stim_labels <- c(KP = "Kleb.", PP = "Pneumo.", LPS = "LPS")
marker_pretty <- c(
  "IL_1beta" = "IL-1β", "IL_1RA"="IL-1RA", "IL_8_1"="IL-8",
  "IL_6"="IL-6", "IL_10"="IL-10", "TNF"="TNF",
  "CCL2"="CCL2", "CCL3"="CCL3", "CCL4"="CCL4"
)

effsize_df <- effsize_df %>%
  mutate(
    stimulation = fct_relabel(stimulation, ~ stim_labels[.x] %||% .x),
    marker = factor(marker_pretty[as.character(marker)] %||% as.character(marker))
  )

# 全局排序
marker_order   <- effsize_df %>% arrange(desc(abs(g))) %>% pull(marker) %>% unique()
contrast_order <- c("CTS1 vs CTS2","CTS1 vs CTS3","CTS2 vs CTS3")

effsize_df <- effsize_df %>%
  mutate(marker = factor(marker, levels = marker_order),
         contrast_lab = factor(contrast_lab, levels = contrast_order))

## ===== 条纹背景 =====
strip_bg <- effsize_df %>%
  group_by(stimulation, contrast_lab) %>%
  distinct(marker) %>%
  arrange(stimulation, contrast_lab, marker) %>%
  mutate(row_id = row_number(),
         ymin = row_id - 0.5, ymax = row_id + 0.5) %>%
  filter(row_id %% 2 == 1) %>%
  ungroup()

## ===== 作图（红=第一组更高，蓝=第二组更高） =====
color_mapping <- c(
  "ns"="#9a9a9a",
  "* UP"="#ffb3b3","** UP"="#ff8080","*** UP"="#ff4d4d","**** UP"="#e60000",
  "* DOWN"="#b3c6ff","** DOWN"="#809fff","*** DOWN"="#4d79ff","**** DOWN"="#0040ff"
)

p_all <- ggplot(effsize_df, aes(x = g, y = marker)) +
  geom_rect(data = strip_bg, inherit.aes = FALSE,
            aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
            fill = "gray95", colour = NA) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.4) +
  geom_segment(aes(x = conf.low, xend = conf.high, yend = marker), size = 0.45) +
  geom_point(aes(fill = signif_dir), shape = 22, size = 3) +
  geom_text(
    data = effsize_df %>% group_by(stimulation, contrast_lab) %>%
      mutate(x_lab = max(conf.high, na.rm = TRUE) + 0.10) %>% ungroup(),
    aes(x = x_lab, label = ifelse(p_BH.signif == "ns","",p_BH.signif)),
    hjust = 0, size = 3.2, family = "sans"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.35))) +
  scale_y_discrete(expand = expansion(add = c(0.5, 0.5))) +
  scale_fill_manual(values = color_mapping, na.value = "#9a9a9a") +
  facet_grid(rows = vars(stimulation), cols = vars(contrast_lab), scales = "free_x") +
  labs(title = "Pairwise comparisons of cytokine responses between CTS subtypes",
       x = "Hedges’ g (positive = higher in first group)", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y   = element_text(face = "bold"),
    strip.text.x  = element_text(face = "bold"),
    strip.text.y  = element_text(face = "bold"),
    panel.grid    = element_blank(),
    axis.line.x   = element_line(size = 0.25),
    legend.position = "none",
    plot.margin   = margin(10, 30, 10, 10),
    panel.spacing.x = unit(10, "pt"),
    panel.spacing.y = unit(8, "pt")
  )

print(p_all)
ggsave("Figure/CTS_effectsize_first_minus_second.svg", p_all, width = 12, height = 8, dpi = 300)

## =====（可选）一致性自检：方向是否与均值差一致？ =====
check_dir <- effsize_df %>%
  group_by(stimulation, contrast_lab, marker) %>%
  summarise(sign_g = sign(first(g)),
            .groups = "drop")
# 若全部在 {-1,0,1}，即可；需要进一步严检可对照前面 md$delta
