##
Count_data <- import("L:/basic/divg/CEMM-Infectious disease/Erik Michels/Alle data/Available data/Gene expression EB&Optimact/OPTIMACT_EB genes.csv")
Patients <- import("L:/basic/divg/CEMM-Infectious disease/Erik Michels/Alle data/Available data/Gene expression EB&Optimact/SK_Amsterdam.csv")
sampleInfo <- select(Patients, "Sample", "RIN", "Sample ID","RNAseq")
sampleInfo <- filter(sampleInfo, sampleInfo$RNAseq == 1)

##
sampleInfo <- sampleInfo[!is.na(sampleInfo$RIN), ]
sampleInfo <- filter(sampleInfo, sampleInfo$RIN >= 5)

row.names(Count_data) <- Count_data$V1
Count_data$V1 <- NULL
Count_data <- Count_data[ , colnames(Count_data) %in% sampleInfo$Sample]


##
healthies <- import("L:/basic/divg/CEMM-Infectious disease/Erik Michels/Alle data/Available data/Gene expression EB&Optimact/Combined_counts_QNS31079_HConly.csv")
row.names(healthies) <- healthies$V1
healthies$V1 <- NULL
colnames(healthies) <- paste(colnames(healthies), ".HV")
colnames(healthies) <- gsub(" ", "", colnames(healthies))
Count_data <- merge(Count_data, healthies, by = 0)
row.names(Count_data) <- Count_data$Row.names
Count_data$Row.names <- NULL

##
sampleInfo2 <- import("L:/basic/divg/CEMM-Infectious disease/Erik Michels/Alle data/Available data/Gene expression EB&Optimact/QNS31079_HCs.csv")
sampleInfo2 <- select(sampleInfo2, "Sample", "RIN", "Sample ID")
sampleInfo2$Sample <- paste(sampleInfo2$Sample, ".HV")
sampleInfo2$Sample <- gsub(" ", "", sampleInfo2$Sample) 
sampleInfo$RNAseq <- NULL
sampleInfo <- rbind(sampleInfo2, sampleInfo)

sampleInfo <- sampleInfo[order(match(sampleInfo$Sample, colnames(Count_data))), ]
sum(sampleInfo$Sample == colnames(Count_data))

# Load required libraries
p_load(DESeq2)

# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(Count_data, colData = sampleInfo, design = ~1)

# Perform size factor estimation and normalization
dds <- DESeq(dds)

# Get normalized counts
normalized_counts <- counts(dds, normalized = TRUE)

# Filtering lowly expressed genes 
keep <- rowSums(counts(dds, normalized = TRUE) > 10) >= 3
dds <- dds[keep,]

# Quality control and visualization (you can customize this as needed)
vsd <- varianceStabilizingTransformation(dds)

# Identify and label the outliers
# Create a PCA plot with labels for RIN outliers
pca_result <- plotPCA(vsd, intgroup = "RIN")

## check stange patients
vsd$group <- ifelse(grepl("EB", vsd$`Sample ID`), "EB", 
                    ifelse(grepl("ELDER", vsd$`Sample ID`), "EB", 
                           ifelse(grepl("BO", vsd$`Sample ID`), "OPT", 
                                  ifelse(grepl("PANAMO|panamo", vsd$`Sample ID`), "COV_ICU", 
                                         ifelse(grepl("covid|COVID", vsd$`Sample ID`), "COV_BIO", "healthy")))))
vsd$group <- ifelse(vsd$`Sample ID` == 1190, "EB", vsd$group)
vsd$group <- ifelse(vsd$`Sample ID` == 3141, "EB", vsd$group)
vsd$group <- ifelse(vsd$`Sample ID` == 2080, "OPT", vsd$group)

##
table(vsd$group)

##
plotPCA(vsd, intgroup = "group")


##
##sample <- vsd$`Sample ID`
##group<- vsd$group
##both <- cbind(sample, group) %>% as.data.frame()
##both <- filter(both, both$group == "pizza")

##check <- vsd[vsd$group == "pizza", ] %>%
##check$`Sample ID`
##head(check)
##

##
vsd$`Sample ID` <- gsub("EB_3.0 NL57847.018.16", "", vsd$`Sample ID`)
vsd$`Sample ID` <- gsub("EB_2.0 NL57847.018.16", "", vsd$`Sample ID`)
vsd$`Sample ID` <- gsub("ELDER-BIOME_2.0 NL57847.018.16", "", vsd$`Sample ID`)
vsd$`Sample ID` <- gsub("ELDER-BIOME_3.0 NL57847.018.16", "", vsd$`Sample ID`)
vsd$`Sample ID` <- gsub("EB NL57847.018.16", "", vsd$`Sample ID`)
vsd$`Sample ID` <- gsub("_d0", "", vsd$`Sample ID`)
vsd$`Sample ID` <- gsub("_D0", "", vsd$`Sample ID`)
vsd$`Sample ID` <- gsub("PAXGENE", "", vsd$`Sample ID`)
vsd$`Sample ID` <- gsub("_DAY 0", "", vsd$`Sample ID`)
vsd$`Sample ID` <- gsub("_DAY0", "", vsd$`Sample ID`)
vsd$`Sample ID` <- gsub(" ", "", vsd$`Sample ID`)
vsd$`Sample ID` <- gsub("EBNL57847.018.21", "", vsd$`Sample ID`)
vsd$`Sample ID` <- gsub("BONL57923.018.16", "", vsd$`Sample ID`)
vsd$`Sample ID` <- gsub("SG", "", vsd$`Sample ID`)
vsd$`Sample ID` <- ifelse(vsd$group == "OPT", paste(vsd$`Sample ID`, "OPT"), vsd$`Sample ID`)

##
setwd("L:/basic/divg/CEMM-Infectious disease/Erik Michels/Alle data/Available data/Gene expression EB&Optimact")
vsd_RNAseq_CAP_ward <- vsd
save(vsd_RNAseq_CAP_ward, file = "vsd_RNAseq_CAP_ward.Rdata")

##
vsd_RNAseq_CAP_ward$`Sample ID`

## get ID evelien
evelien <- import("C:/Users/P074130/OneDrive - Amsterdam UMC/Bureaublad/ID.xlsx") %>% as.data.frame()
vsd_RNAseq_CAP_ward$`Sample ID` <- paste("EB_", vsd_RNAseq_CAP_ward$`Sample ID`)
vsd_RNAseq_CAP_ward$`Sample ID` <- gsub(" ", "", vsd_RNAseq_CAP_ward$`Sample ID`)
vsd_RNAseq_CAP_ward <- vsd_RNAseq_CAP_ward[vsd_RNAseq_CAP_ward$`Sample ID` %in% evelien]


##
##
p_load(SepstratifieRSRS <- assay(vsd)
rownames(SRS) <- sub("\\..*$", "", rownames(SRS))
predictions <- stratifyPatients(t(SRS))
predictions_extended <- stratifyPatients(t(SRS), gene_set = "extended")

plotAlignedSamples(predictions)
plotAlignedSamples(predictions_extended)

short_SRSq <- predictions@SRSq
long_SRSq <- predictions_extended@SRSq

table(predictions@SRS)
table(predictions_extended@SRS)

short_SRS_class <- predictions@SRS %>% as.data.frame()
long_SRS_class <- predictions_extended@SRS %>% as.data.frame()
SRS_scores <- merge(short_SRSq, long_SRSq, by = 0)
SRS_scores <- merge(SRS_scores, short_SRS_class, by.x = "Row.names", by.y = 0)
SRS_scores <- merge(SRS_scores, long_SRS_class, by.x = "Row.names", by.y = 0)
colnames(SRS_scores) <- c("Row.names", "SRSq_short", "SRSq_long", "SRS_7", "SRS_19")

##


##
name <- vsd$Sample
name2 <- vsd$`Sample ID`
identifier <- cbind(name, name2) %>% as.data.frame()
SRS_scores <- merge(identifier, SRS_scores, by.x = "name", by.y = "Row.names")




