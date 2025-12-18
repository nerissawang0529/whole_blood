rm(list = ls())

## install.packages("pacman")
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
               tableone, kmed, survival, survminer, skimr, lattice, reshape2, grDevices, ggrepel)
install.packages(c("cluster", "ConsensusClusterPlus", "factoextra", "ggplot2"))
remove.packages("cli", lib=.libPaths()[1])
install.packages("cli", dependencies=TRUE)
install.packages("ComplexHeatmap")



if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install ConsensusClusterPlus from Bioconductor
BiocManager::install("ConsensusClusterPlus")

library(ConsensusClusterPlus)
library(cluster)  # For Gower distance
library(ConsensusClusterPlus)  # For Consensus Clustering
library(factoextra)  # For visualization
library(ggplot2)

library(dplyr)
library(MASS)
library(magrittr)
library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)
library(readxl)

studydata_patients <- read.csv("Original_data/studydata_patients.csv")
studydata_HV <- read.csv("Original_data/studydata_HV.csv")
studydata_combined <- rbind(studydata_patients, studydata_HV)

Luminex_CAP_COVID_HV <- read.csv("Original_data/Luminex_CAP_COVID_HV.csv")
merged_data <- merge(Luminex_CAP_COVID_HV, studydata_patients, 
                     by.x = "ID", by.y = "EB_id")

##export form
destination_folder <- "Original_data/" 
export_file_name <- "whole_blood_stimulation.csv" 
write.csv(merged_data, file.path(destination_folder, export_file_name), row.names = FALSE)
cat("File exported to:", file.path(destination_folder, export_file_name), "\n")



#merged_data <- merged_data %>% filter(stimulation != "M")

# Convert categorical variables into dummy variables
merged_data <- merged_data %>%
  mutate(across(c("stimulation"), as.factor))

# Convert categorical variables to dummy variables
dummy_vars <- model.matrix(~ stimulation - 1, data = merged_data) %>%
  as.data.frame()

# Merge dummy variables back into the original dataset
merged_data <- cbind(merged_data %>% dplyr::select(-stimulation), dummy_vars)

# Ensure all cytokine markers are numeric
#num_markers <- c("minis_M_CCL2", "minis_M_CCL4", "minis_M_IL_1RA", "minis_M_IL_8_1", "minis_M_TNF", "minis_M_CCL3", "minis_M_IL_1beta", "minis_M_IL_6", "minis_M_IL_10")
num_markers <- c("CCL2","CCL4","IL_1RA","IL_8_1","TNF","CCL3","IL_1beta","IL_6","IL_10")

merged_data[num_markers] <- lapply(merged_data[num_markers], function(x) as.numeric(as.character(x)))

# **Scale the numeric cytokine markers**
# or log the data
#merged_data[num_markers] <- scale(merged_data[num_markers])
merged_data[, num_markers] <- log10(merged_data[, num_markers])
# Ensure binary categorical variables are treated as factors
cat_vars <- c("stimulationKP", "stimulationLPS", "stimulationPP", "stimulationM")

merged_data[cat_vars] <- lapply(merged_data[cat_vars], as.factor)

clustering_data <- merged_data[, c(num_markers, cat_vars)]

# Compute Gower distance (handles mixed data types)
distance_matrix <- daisy(clustering_data, metric = "gower")

# Run Consensus Clustering
results <- ConsensusClusterPlus(
  as.matrix(distance_matrix),  # Convert distance matrix to numeric format
  maxK = 6,  # Test clusters from k = 2 to 6
  reps = 100,  # Number of resampling iterations
  pItem = 0.8,  # 80% of samples used in each iteration
  pFeature = 1.0,  # Use all features
  clusterAlg = "hc",  # Hierarchical clustering
  distance = "euclidean",  # Required by ConsensusClusterPlus
  seed = 1234,  # Set seed for reproducibility
  title = "ConsensusClustering"
)

# Select the best number of clusters (for example, k = 3)
best_k <- 3  # Change if needed based on visualizations
clusters <- results[[best_k]]$consensusClass
merged_data$Cluster <- as.factor(clusters)
table(merged_data$Cluster)















merged_data_clean <- na.omit(merged_data[, num_markers])
merged_data <- merged_data[complete.cases(merged_data[, num_markers]), ]

pca_res <- prcomp(merged_data_clean, center = TRUE, scale. = TRUE)

autoplot(pca_res, data = merged_data, colour = "Cluster") +
  ggtitle("PCA Plot Colored by Consensus Clustering") +
  theme_minimal()




















# Visualize the consensus matrix heatmap
heatmap(results[[best_k]]$consensusMatrix, main = "Consensus Clustering Heatmap")

# Plot cumulative distribution function (CDF)
# Extract cumulative distribution function (CDF) values from results
k_values <- 2:length(results)  # Extract tested cluster numbers
cdf_values <- lapply(k_values, function(k) results[[k]]$consensusMatrix)

# Convert CDF values into a data frame
cdf_df <- data.frame(
  k = rep(k_values, each = length(cdf_values[[1]])),
  consensus_index = unlist(lapply(cdf_values, function(mat) as.vector(mat)))
)
# Load ggplot2 if not already loaded
library(ggplot2)

# Create the CDF plot
ggplot(cdf_df, aes(x = consensus_index, color = as.factor(k))) +
  stat_ecdf(geom = "step", size = 1) +
  labs(title = "Consensus CDF Plot", x = "Consensus Index", y = "CDF",
       color = "Number of Clusters (k)") +
  theme_minimal()


# Load required library
library(factoextra)
# Perform hierarchical clustering on Gower distance
hc <- hclust(distance_matrix, method = "average")
# Plot dendrogram with clusters
fviz_dend(hc, k = best_k, rect = TRUE, rect_fill = TRUE, rect_border = "jco", cex = 0.6)

# Assign cluster labels from ConsensusClusterPlus results
merged_data$Cluster <- as.factor(results[[best_k]]$consensusClass)

# Boxplot example: IL-6 levels across clusters
ggplot(merged_data, aes(x = Cluster, y = minis_M_IL_6, fill = Cluster)) +
  geom_boxplot() +
  labs(title = "IL-6 Levels by Cluster", x = "Cluster", y = "IL-6 (Z-score)") +
  theme_minimal()





# Sample only 100 rows to make the plot readable
set.seed(123)
sample_indices <- sample(1:nrow(consensus_matrix), 100)
sub_matrix <- consensus_matrix[sample_indices, sample_indices]

hc <- hclust(dist(sub_matrix), method = "ward.D2")
plot(hc, main = "Hierarchical Clustering (Subset)")

