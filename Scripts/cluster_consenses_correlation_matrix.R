### make matrix of piD versus "markers" (cytokine_stimulation combinations n=36)
rm(list = ls())

studydata_patients <- read.csv("Original_data/studydata_patients.csv")
#studydata_HV <- read.csv("Original_data/studydata_HV.csv")
#studydata_combined <- rbind(studydata_patients, studydata_HV)

Luminex_CAP_COVID_HV <- read.csv("Original_data/Luminex_CAP_COVID_HV.csv")
merged_data <- merge(Luminex_CAP_COVID_HV, studydata_patients, 
                     by.x = "ID", by.y = "EB_id")

mdata <- merged_data[,1:12]
mdata$Day <- NULL

library(reshape2)
ldata <- reshape2::melt(mdata, id=c("ID","stimulation") )
ldata$level <- log10(ldata$value)

ldata$value <- NULL


raw_w <- reshape2::dcast(ldata, ID ~ stimulation + variable, value.var="level" )

row.names(raw_w) <- raw_w$ID
raw_w$ID <- NULL



raw_w_t <- t(raw_w)
raw_w_t <- as.data.frame(raw_w_t)
colnames(raw_w_t) <- raw_w_t[1, ]  # Set first row as column names
raw_w_t <- raw_w_t[-1, ]           # Remove the first row
clustering_data <- raw_w_t
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
