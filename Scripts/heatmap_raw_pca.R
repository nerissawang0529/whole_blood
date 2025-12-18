library(pheatmap)
library(RColorBrewer)

# === 1. Combine PC2 and selected biomarkers ===
heatmap_data <- data.frame(
  PC2           = pca_df$PC2,
  TNF_LPS       = ldata_wide2$TNF_LPS,
  TNF_KP        = ldata_wide2$TNF_KP,
  TNF_PP        = ldata_wide2$TNF_PP,
  IL_1RA_LPS    = ldata_wide2$IL_1RA_LPS,
  IL_1RA_KP     = ldata_wide2$IL_1RA_KP,
  IL_1RA_PP     = ldata_wide2$IL_1RA_PP,
  IL_1beta_LPS  = ldata_wide2$IL_1beta_LPS,
  IL_1beta_KP   = ldata_wide2$IL_1beta_KP,
  IL_1beta_PP   = ldata_wide2$IL_1beta_PP,
  IL_6_LPS      = ldata_wide2$IL_6_LPS,
  IL_6_KP       = ldata_wide2$IL_6_KP,
  IL_6_PP       = ldata_wide2$IL_6_PP
)

# === 2. Order patients by PC2 ===
heatmap_data_ordered <- heatmap_data[order(heatmap_data$PC2), ]

# === 3. Extract and scale biomarker matrix (Z-score normalization) ===
biomarker_matrix <- as.matrix(heatmap_data_ordered[, -1])  # remove PC2
biomarker_matrix <- scale(biomarker_matrix)

# === 4. Set custom marker order (not based on correlation) ===
marker_order <- c(
  "TNF_LPS", "TNF_KP", "TNF_PP",
  "IL_1RA_LPS", "IL_1RA_KP", "IL_1RA_PP",
  "IL_1beta_LPS", "IL_1beta_KP", "IL_1beta_PP",
  "IL_6_LPS", "IL_6_KP", "IL_6_PP"
)
toplot <- t(biomarker_matrix)[marker_order, ]

# === 5. Create PC2 annotation bar on top ===
annotation_col <- data.frame(PC2 = heatmap_data_ordered$PC2)
rownames(annotation_col) <- rownames(heatmap_data_ordered)

pc2_colors <- colorRampPalette(c("blue", "white", "red"))(100)
annotation_colors <- list(PC2 = pc2_colors)

# === 6. Create symmetric color scale ===
custom_palette <- colorRampPalette(c("blue", "grey90", "red"))(200)
max_abs <- max(abs(biomarker_matrix), na.rm = TRUE)
breaks <- seq(-max_abs, max_abs, length.out = 201)

# === 7. Draw the heatmap ===
pheatmap(
  toplot,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  color = custom_palette,
  breaks = breaks,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  main = "Cytokine Responses (LPS, Kneu, Pneu) Ordered by PC2"
)





library(pheatmap)
library(RColorBrewer)

# === 1. Combine PC1 and selected biomarkers ===
heatmap_data <- data.frame(
  PC1           = pca_df$PC1,
  TNF_LPS       = ldata_wide2$TNF_LPS,
  TNF_KP        = ldata_wide2$TNF_KP,
  TNF_PP        = ldata_wide2$TNF_PP,
  IL_1RA_LPS    = ldata_wide2$IL_1RA_LPS,
  IL_1RA_KP     = ldata_wide2$IL_1RA_KP,
  IL_1RA_PP     = ldata_wide2$IL_1RA_PP,
  IL_1beta_LPS  = ldata_wide2$IL_1beta_LPS,
  IL_1beta_KP   = ldata_wide2$IL_1beta_KP,
  IL_1beta_PP   = ldata_wide2$IL_1beta_PP,
  IL_6_LPS      = ldata_wide2$IL_6_LPS,
  IL_6_KP       = ldata_wide2$IL_6_KP,
  IL_6_PP       = ldata_wide2$IL_6_PP
)

# === 2. Order patients by PC1 ===
heatmap_data_ordered <- heatmap_data[order(heatmap_data$PC1), ]

# === 3. Extract and scale biomarker matrix (Z-score normalization) ===
biomarker_matrix <- as.matrix(heatmap_data_ordered[, -1])  # remove PC1
biomarker_matrix <- scale(biomarker_matrix)

# === 4. Set custom marker order (not based on correlation) ===
marker_order <- c(
  "TNF_LPS", "TNF_KP", "TNF_PP",
  "IL_1RA_LPS", "IL_1RA_KP", "IL_1RA_PP",
  "IL_1beta_LPS", "IL_1beta_KP", "IL_1beta_PP",
  "IL_6_LPS", "IL_6_KP", "IL_6_PP"
)
toplot <- t(biomarker_matrix)[marker_order, ]

# === 5. Create PC1 annotation bar on top ===
annotation_col <- data.frame(PC1 = heatmap_data_ordered$PC1)
rownames(annotation_col) <- rownames(heatmap_data_ordered)

pc1_colors <- colorRampPalette(c("blue", "white", "red"))(100)
annotation_colors <- list(PC1 = pc1_colors)

# === 6. Create symmetric color scale ===
custom_palette <- colorRampPalette(c("blue", "grey90", "red"))(200)
max_abs <- max(abs(biomarker_matrix), na.rm = TRUE)
breaks <- seq(-max_abs, max_abs, length.out = 201)

# === 7. Draw the heatmap ===
pheatmap(
  toplot,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  color = custom_palette,
  breaks = breaks,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  main = "Cytokine Responses (LPS, Kneu, Pneu) Ordered by PC1"
)






library(pheatmap)
library(RColorBrewer)

# === 1. Create the full data frame ===
heatmap_data <- data.frame(
  TNF_LPS       = ldata_wide2$TNF_LPS,  # Used only for ordering & annotation
  TNF_KP        = ldata_wide2$TNF_KP,
  TNF_PP        = ldata_wide2$TNF_PP,
  IL_1RA_LPS    = ldata_wide2$IL_1RA_LPS,
  IL_1RA_KP     = ldata_wide2$IL_1RA_KP,
  IL_1RA_PP     = ldata_wide2$IL_1RA_PP,
  IL_1beta_LPS  = ldata_wide2$IL_1beta_LPS,
  IL_1beta_KP   = ldata_wide2$IL_1beta_KP,
  IL_1beta_PP   = ldata_wide2$IL_1beta_PP,
  IL_6_LPS      = ldata_wide2$IL_6_LPS,
  IL_6_KP       = ldata_wide2$IL_6_KP,
  IL_6_PP       = ldata_wide2$IL_6_PP
)

# === 2. Order patients by raw TNF_LPS ===
heatmap_data_ordered <- heatmap_data[order(heatmap_data$TNF_LPS), ]

# === 3. Extract biomarker matrix WITHOUT TNF_LPS and apply Z-score ===
biomarker_matrix <- as.matrix(heatmap_data_ordered[, -1])  # Remove TNF_LPS
biomarker_matrix <- scale(biomarker_matrix)

# === 4. Create annotation bar with raw TNF_LPS ===
annotation_col <- data.frame(TNF_LPS = heatmap_data_ordered$TNF_LPS)
rownames(annotation_col) <- rownames(heatmap_data_ordered)

tnf_colors <- colorRampPalette(c("blue", "white", "red"))(100)
annotation_colors <- list(TNF_LPS = tnf_colors)

# === 5. Symmetric color palette for scaled data ===
custom_palette <- colorRampPalette(c("blue", "grey90", "red"))(200)
max_abs <- max(abs(biomarker_matrix), na.rm = TRUE)
breaks <- seq(-max_abs, max_abs, length.out = 201)

# === 6. Draw the heatmap ===
pheatmap(
  t(biomarker_matrix),               # Markers as rows, patients as columns
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  color = custom_palette,
  breaks = breaks,
  main = "Cytokine Responses Ordered by Raw TNF_LPS"
)

