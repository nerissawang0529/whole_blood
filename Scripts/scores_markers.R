DIP_continous <- read.csv("Original_data/DIP_predictions_20250606_continous.csv")
DIP_continous$TREM_1 <- NULL
DIP_continous$IL_6 <- NULL
DIP_continous$Procalcitonin <- NULL
studydata_patients <- read.csv("Original_data/studydata_patients.csv")
Luminex_CAP_COVID_HV <- read.csv("Original_data/Luminex_CAP_COVID_HV.csv")
merged_data <- merge(Luminex_CAP_COVID_HV, studydata_patients, 
                     by.x = "ID", by.y = "EB_id")
combined_data <- merge(DIP_continous, merged_data, by = "ID")
selected_data <- combined_data[, c("ID", "cDIP", "CCL2", "CCL4", "IL_1RA", "IL_8_1", "TNF", "CCL3", "IL_1beta",
                                   "IL_6", "IL_10", "stimulation", "group", "MEWS_score", "CURB_score",
                                   "PSI_new", "mortality_d30", "mortality_d90")]

library(dplyr)
library(tidyr)

library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)

# Define the markers you want to include
markers <- c("TNF", "IL_1RA", "IL_1beta", "IL_6", "IL_8_1", "IL_10", "CCL2", "CCL3", "CCL4")

# Reshape selected_data to wide format: marker_stim rows, patient IDs as columns
data_wide <- selected_data %>%
  filter(stimulation != "M") %>%  # 🔥 Exclude "_M"
  select(ID, stimulation, all_of(markers)) %>%
  pivot_longer(cols = all_of(markers), names_to = "marker", values_to = "value") %>%
  mutate(marker_stim = paste(marker, stimulation, sep = "_")) %>%
  select(ID, marker_stim, value) %>%
  pivot_wider(names_from = ID, values_from = value)

# Convert to matrix
heatmap_matrix <- as.data.frame(data_wide)
rownames(heatmap_matrix) <- heatmap_matrix$marker_stim
heatmap_matrix <- heatmap_matrix[, -1]

# Z-score scaling across rows
heatmap_matrix_z <- t(scale(t(as.matrix(heatmap_matrix))))

# === Add DIP_continous annotation ===
# Ensure ID names match those in heatmap columns
dip_order <- DIP_continous %>%
  select(ID, cDIP) %>%
  filter(ID %in% colnames(heatmap_matrix_z)) %>%
  arrange(cDIP)

# Reorder heatmap columns by cDIP
heatmap_matrix_z <- heatmap_matrix_z[, as.character(dip_order$ID)]

# Annotation bar
annotation_col <- data.frame(cDIP = dip_order$cDIP)
rownames(annotation_col) <- as.character(dip_order$ID)

# Color palette for annotation bar
cDIP_colors <- colorRampPalette(c("white", "darkred"))(100)
annotation_colors <- list(cDIP = cDIP_colors)

# Color palette and breaks for heatmap
custom_palette <- colorRampPalette(c("blue", "grey90", "red"))(200)
max_abs <- max(abs(heatmap_matrix_z), na.rm = TRUE)
breaks <- seq(-3, 3, length.out = 201)



# Extract marker names from rownames (e.g., "TNF_LPS" → "TNF")
marker_base_order <- gsub("_(LPS|KP|PP|M)$", "", rownames(heatmap_matrix_z))

# Create a factor that forces the desired grouping order
marker_stimulation_order <- order(factor(marker_base_order, levels = markers))

# Reorder heatmap matrix rows
heatmap_matrix_z <- heatmap_matrix_z[marker_stimulation_order, ]

# === Draw the heatmap ===
pheatmap(
  heatmap_matrix_z,
  cluster_rows = FALSE,         # No row clustering
  cluster_cols = FALSE,         # No column clustering
  show_colnames = FALSE,        # Hide column labels (patient IDs)
  show_rownames = TRUE,         # Show marker_stim labels
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  color = custom_palette,
  breaks = breaks,
  main = "Marker-Stimulation Cytokine Responses Ordered by cDIP"
)



library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)

# === Step 1: Define markers ===
markers <- c("TNF", "IL_1RA", "IL_1beta", "IL_6", "IL_8_1", "IL_10", "CCL2", "CCL3", "CCL4")

# === Step 2: Reshape to wide format with ID_stimulation as columns ===
data_long <- selected_data %>%
  filter(stimulation != "M") %>%
  pivot_longer(cols = all_of(markers), names_to = "marker", values_to = "value") %>%
  mutate(
    marker_stim = paste(marker, stimulation, sep = "_"),
    ID_stimulation = paste(ID, stimulation, sep = "_")
  )

data_wide <- data_long %>%
  select(ID_stimulation, marker_stim, value) %>%
  pivot_wider(names_from = ID_stimulation, values_from = value)

# === Step 3: Convert to matrix & scale ===
heatmap_matrix <- as.data.frame(data_wide)
rownames(heatmap_matrix) <- heatmap_matrix$marker_stim
heatmap_matrix <- heatmap_matrix[, -1]
heatmap_matrix_z <- t(scale(t(as.matrix(heatmap_matrix))))

# === Step 4: Create annotation & sort columns by MEWS_score ===
annotation_col_mews <- selected_data %>%
  filter(stimulation != "M") %>%
  mutate(ID_stimulation = paste(ID, stimulation, sep = "_")) %>%
  filter(ID_stimulation %in% colnames(heatmap_matrix_z)) %>%
  distinct(ID_stimulation, .keep_all = TRUE) %>%
  arrange(MEWS_score)  # ✅ sort by MEWS_score from low to high

# Reorder heatmap columns
heatmap_matrix_z <- heatmap_matrix_z[, annotation_col_mews$ID_stimulation]

# Create annotation data frame
annotation_col <- data.frame(MEWS_score = annotation_col_mews$MEWS_score)
rownames(annotation_col) <- annotation_col_mews$ID_stimulation

# === Step 5: Marker row ordering ===
marker_base_order <- gsub("_(LPS|KP|PP)$", "", rownames(heatmap_matrix_z))
marker_stimulation_order <- order(factor(marker_base_order, levels = markers))
heatmap_matrix_z <- heatmap_matrix_z[marker_stimulation_order, ]

# === Step 6: Draw the heatmap ===
custom_palette <- colorRampPalette(c("blue", "grey90", "red"))(200)
breaks <- seq(-3, 3, length.out = 201)

pheatmap(
  heatmap_matrix_z,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  show_rownames = TRUE,
  annotation_col = annotation_col,
  color = custom_palette,
  breaks = breaks,
  main = "Cytokine Responses by ID_stimulation (Sorted by MEWS Score)"
)
