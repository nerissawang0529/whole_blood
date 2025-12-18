rm(list = ls())

# Load required packages
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
               tableone, kmed, survival, survminer, skimr, lattice, reshape2, grDevices, ggrepel,
               cluster, ConsensusClusterPlus, factoextra, ComplexHeatmap)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ConsensusClusterPlus")

# Load Data
studydata_patients <- read.csv("Original_data/studydata_patients.csv")
studydata_HV <- read.csv("Original_data/studydata_HV.csv")
studydata_combined <- rbind(studydata_patients, studydata_HV)

Luminex_CAP_COVID_HV <- read.csv("Original_data/Luminex_CAP_COVID_HV.csv")
merged_data <- merge(Luminex_CAP_COVID_HV, studydata_combined, 
                     by.x = "ID", by.y = "EB_id")

# Export merged data
destination_folder <- "Original_data/"
export_file_name <- "whole_blood_stimulation.csv"
write.csv(merged_data, file.path(destination_folder, export_file_name), row.names = FALSE)
cat("File exported to:", file.path(destination_folder, export_file_name), "\n")

# Convert categorical variables to factors
merged_data <- merged_data %>%
  mutate(across(c("stimulation"), as.factor))

# Create dummy variables for stimulation
dummy_vars <- model.matrix(~ stimulation - 1, data = merged_data) %>%
  as.data.frame()
merged_data <- cbind(merged_data %>% dplyr::select(-stimulation), dummy_vars)

# Numeric cytokine markers
num_markers <- c("CCL2","CCL4","IL_1RA","IL_8_1","TNF","CCL3","IL_1beta","IL_6","IL_10")

# Convert to numeric & log-transform cytokine markers
merged_data[num_markers] <- lapply(merged_data[num_markers], function(x) log10(as.numeric(as.character(x))))

# Scale cytokines for clustering
merged_data[num_markers] <- scale(merged_data[num_markers])

# Categorical binary stimulation variables
cat_vars <- c("stimulationKP", "stimulationLPS", "stimulationPP", "stimulationM")
merged_data[cat_vars] <- lapply(merged_data[cat_vars], as.factor)

# Select data for clustering
clustering_data <- merged_data[, c(num_markers, cat_vars)]

# Compute Spearman Correlation-Based Distance
cor_matrix <- cor(t(clustering_data[, num_markers]), method = "spearman")  # Patients vs. Cytokine Response
distance_matrix <- as.dist(1 - cor_matrix)  # Convert to distance

# Run Consensus Clustering
results <- ConsensusClusterPlus(
  as.matrix(distance_matrix),  # Convert distance matrix to numeric format
  maxK = 6,  # Test clusters from k = 2 to 6
  reps = 100,  # Resampling iterations
  pItem = 0.8,  # 80% of samples used in each iteration
  pFeature = 1.0,  # Use all features
  clusterAlg = "hc",  # Hierarchical clustering
  distance = "euclidean",  # Required by ConsensusClusterPlus
  seed = 1234,  # For reproducibility
  title = "ConsensusClustering"
)

# Select the best number of clusters (e.g., k = 3 or k = 5)
best_k <- 3  # Change based on visualization
clusters <- results[[best_k]]$consensusClass
merged_data$Cluster <- as.factor(clusters)

# Print cluster distribution
table(merged_data$Cluster)

# Plot clustering result
fviz_cluster(list(data = clustering_data[, num_markers], cluster = merged_data$Cluster))
