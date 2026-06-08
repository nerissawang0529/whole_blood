############################################################
## Clean workspace
############################################################
rm(list = ls())

############################################################
## Packages
############################################################
library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(svglite)

############################################################
## Load data
############################################################
studydata <- read_csv(
  "Original_data/studydata_with_ms1_il1r2.csv",
  show_col_types = FALSE
)

############################################################
## Global style (LOCKED)
############################################################
FONT_FAMILY <- "Arial"  # Font family used for all text
BASE_SIZE   <- 14       # Base font size (pt)
TITLE_SIZE  <- 14       # Plot title size (pt)
TAG_SIZE    <- 16       # Panel tag size (pt) - NOT bold
AXIS_TITLE  <- 13       # Axis title size (pt)
AXIS_TEXT   <- 11       # Axis tick label size (pt)

LINE_AXIS   <- 0.6      # Axis line width
LINE_TICKS  <- 0.6      # Tick line width
BIN_EDGE_W  <- 0.2      # Histogram bin edge width
MARGINS_MM  <- c(5.5, 5.5, 5.5, 5.5)  # Plot margins (t, r, b, l)

############################################################
## Colors (LOCKED)
############################################################
col_ms1   <- "#2F4B7C"
col_il1r2 <- "#A23B3B"

############################################################
## Theme (LOCKED)
############################################################
base_theme <- theme_classic(base_size = BASE_SIZE) +
  theme(
    text        = element_text(family = FONT_FAMILY),
    plot.title  = element_text(
      size = TITLE_SIZE,
      face = "plain",
      hjust = 0
    ),
    plot.tag    = element_text(
      size = TAG_SIZE,
      face = "plain"
    ),
    axis.title  = element_text(size = AXIS_TITLE),
    axis.text   = element_text(size = AXIS_TEXT),
    axis.line   = element_line(linewidth = LINE_AXIS),
    axis.ticks  = element_line(linewidth = LINE_TICKS),
    plot.margin = margin(
      MARGINS_MM[1],
      MARGINS_MM[2],
      MARGINS_MM[3],
      MARGINS_MM[4]
    )
  )

############################################################
## Histogram settings (LOCKED)
############################################################
n_bins <- 12

############################################################
## MS1 histogram
############################################################
p_ms1 <- ggplot(studydata, aes(x = MS1_percent)) +
  geom_histogram(
    bins      = n_bins,
    fill      = col_ms1,
    color     = "white",
    linewidth = BIN_EDGE_W,
    alpha     = 1,
    na.rm     = TRUE
  ) +
  labs(
    x = "Percentage of MS1 cells",
    y = "Number of patients"
  ) +
  base_theme

############################################################
## IL1R2 histogram
############################################################
p_il1r2 <- ggplot(studydata, aes(x = IL1R2_percent)) +
  geom_histogram(
    bins      = n_bins,
    fill      = col_il1r2,
    color     = "white",
    linewidth = BIN_EDGE_W,
    alpha     = 1,
    na.rm     = TRUE
  ) +
  labs(
    x = "Percentage of IL1R2+ cells",
    y = "Number of patients"
  ) +
  base_theme

############################################################
## Combine panels
## - Locked two-column layout
## - Equal panel widths
## - ERJ-style panel tags: a), b)
## - Use the same layout block for similar two-panel figures
############################################################
final_plot <- p_ms1 + p_il1r2 +
  plot_layout(
    ncol   = 2,
    widths = c(1, 1)
  ) +
  plot_annotation(
    tag_levels = "a",
    tag_suffix = ")"
  )

print(final_plot)

############################################################
## Export as SVG
## - Useful for Illustrator editing
############################################################
ggsave(
  filename = "Figure/Figure_MS1_IL1R2_distribution.svg",
  plot     = final_plot,
  device   = svglite,
  width    = 180,
  height   = 80,
  units    = "mm"
)




#information for il1r2####
x <- studydata$IL1R2_percent

n <- sum(!is.na(x))
zero_n <- sum(x == 0, na.rm = TRUE)
zero_pct <- zero_n / n * 100

med <- median(x, na.rm = TRUE)
q1 <- quantile(x, 0.25, na.rm = TRUE)
q3 <- quantile(x, 0.75, na.rm = TRUE)
min_x <- min(x, na.rm = TRUE)
max_x <- max(x, na.rm = TRUE)

paste0(
  "IL1R2: n = ", n,
  "; zero = ", zero_n, "/", n,
  " (", round(zero_pct, 1), "%)",
  "; median ", round(med, 1), "%",
  " [IQR ", round(q1, 1), "–", round(q3, 1), "]",
  ", range ", round(min_x, 1), "–", round(max_x, 1), "%."
)     




x <- studydata$MS1_percent

n <- sum(!is.na(x))
med <- median(x, na.rm = TRUE)
q1 <- quantile(x, 0.25, na.rm = TRUE)
q3 <- quantile(x, 0.75, na.rm = TRUE)
min_x <- min(x, na.rm = TRUE)
max_x <- max(x, na.rm = TRUE)

paste0(
  "MS1: n = ", n,
  "; median ", round(med, 1), "%",
  " [IQR ", round(q1, 1), "–", round(q3, 1), "]",
  ", range ", round(min_x, 1), "–", round(max_x, 1), "%."
)