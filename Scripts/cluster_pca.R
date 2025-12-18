### make matrix of piD versus "markers" (cytokine_stimulation combinations n=36)


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

raw_w["3005", "PP_IL_1beta"] <- mean(raw_w$PP_IL_1beta, na.rm=T)
raw_w["3350", "PP_CCL3"] <- mean(raw_w$PP_CCL3, na.rm=T)


pca_res <- prcomp(raw_w)

library(ggfortify)

autoplot(pca_res)


autoplot(pca_res, 
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3)


kmeansres <- kmeans(pca_res$x[,1:2], centers = 2) 

pat_clust <- data.frame(cluster=as.factor(kmeansres$cluster))

autoplot(pca_res, data=pat_clust, colour="cluster",
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3)

