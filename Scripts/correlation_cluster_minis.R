rm(list = ls())

studydata_patients <- read.csv("Original_data/studydata_patients.csv")
Luminex_CAP_COVID_HV <- read.csv("Original_data/Luminex_CAP_COVID_HV_log.csv")

merged_data <- merge(Luminex_CAP_COVID_HV, studydata_patients, 
                     by.x = "ID", by.y = "EB_id")

mdata <- merged_data[, c(1, 12, 14:22)]
mdata <- mdata %>% filter (!stimulation == "M")

library(reshape2)
ldata <- reshape2::melt(mdata, id=c("ID","stimulation") )

raw_w <- reshape2::dcast(ldata, ID ~ stimulation + variable, value.var="value" )

row.names(raw_w) <- raw_w$ID
raw_w$ID <- NULL

corm <- cor(raw_w, method = "spearman",
            use = "pairwise.complete.obs")

library(corrplot)
corrplot(corm)
corrplot(corm, method = "color", type = "upper", 
         col = colorRampPalette(c("white", "black"))(300),
         tl.col = "black", tl.srt = 75)


distmat <- as.dist(1 - corm)

plot(distmat)

hco <- hclust(distmat, method = "complete")

plot(hco)

plot(as.dendrogram(hco))

breaks <- seq(-1, 1, length.out = 101) 

library(pheatmap)
pheatmap(corm, 
         clustering_distance_rows = distmat,
         clustering_distance_cols = distmat,
         #Rowv = as.dendrogram(hco), 
         #Colv = as.dendrogram(hco), 
         scale = "none",
         color = colorRampPalette(c("white", "black"))(300),
         breaks = breaks,
         main = "Spearman Correlation Heatmap")

# Cut dendrogram into 3 clusters (or adjust k)
nclust <- cutree(hco, k=5)
table(nclust)

# Convert to a data frame for easy viewing
cluster_assignments <- data.frame(Marker = names(nclust), Cluster = nclust)

# View which markers belong to which cluster
cluster_assignments <- cluster_assignments %>%
  arrange(Cluster)
print(cluster_assignments)
