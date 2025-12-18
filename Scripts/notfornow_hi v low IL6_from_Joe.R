rm(list = ls ())

eb <- read.csv("Original_data/ELDER-BIOME_export_20230605.csv", sep=";")

ebs <- eb[, grepl("outcome_pat_microbiology|Participant.Id|base_date_admission", colnames(eb))]

colnames(ebs)[1] <- "Participant.Id"

library(reshape2)
melted_df <- reshape2::melt(ebs, id.vars = c("Participant.Id", "base_date_admission"))



melted_df$variable <- as.character(melted_df$variable)

melted_df$varx <-  sapply(strsplit(melted_df$variable, "_"), "[", 7)

melted_df$num <-  sapply(strsplit(melted_df$variable, "_"), "[", 6)

melted_df$variable <- NULL

wide_df <- reshape2::dcast(melted_df, Participant.Id + num + base_date_admission ~ varx)

table(wide_df$Sample.origin )


wide_df$base_date_admission <- as.Date(wide_df$base_date_admission, format="%d-%m-%Y")

wide_df$Sample.date <- as.Date(wide_df$Sample.date, format="%d-%m-%Y")


wide_df$adm_samp_dd <- difftime(wide_df$Sample.date, wide_df$base_date_admission, units="days")



bc <- wide_df[wide_df$Sample.origin=="Blood culture",]



wbs <- read.csv("Original_data/whole_blood_stimuli_long.csv")


table(wbs$HiIL6)


dfwbs <- unique(wbs[,c("ID","HiIL6","group")])

table(dfwbs$HiIL6)
table(dfwbs$HiIL6, dfwbs$group)
## patinet key code

key <- read.csv("Original_data/SK_Amsterdam.csv")

count_mat <- read.csv("Original_data/Combined_counts_Amsterdam.csv")

row.names(count_mat) <- sapply(strsplit(count_mat$X, "\\."), "[", 1)

count_mat$X <- NULL



table(colnames(count_mat) %in% key$Sample)

Elder <- key[grepl("EB|ELDER", key$Sample.ID ) & key$RNAseq==1,]


library(stringr)

Elder$ID <- str_match(Elder$Sample.ID, pattern='\\d\\d\\d\\d_')

Elder$ID <- str_remove(Elder$ID , "_")


table(dfwbs$ID %in% Elder$ID)


dfwbseq <- dfwbs[dfwbs$ID %in% Elder$ID,]
table(dfwbseq$group, dfwbseq$HiIL6)




### CAP


dfwbscap <- dfwbseq[dfwbseq$group=="CAP",]

table(dfwbscap$HiIL6)


dfwbscapm <- merge(dfwbscap, Elder, by="ID")

row.names(dfwbscapm) <- dfwbscapm$Sample






#table(bcs_str$BACTERIA.found[bcs_str$Participant.Id %in% CASES])

count_eb <- count_mat[,  dfwbscapm$Sample ]

table(row.names(dfwbscapm) == colnames(count_eb))

# Create a folder to save the packages
dir.create("DESeq2_offline_lib")





dds <- DESeqDataSetFromMatrix(countData = count_eb,
                              colData = dfwbscapm,
                              design= ~ HiIL6)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, contrast=c("HiIL6","HighIL6_M","LowIL6_M"))
# or to shrink log fold changes association with condition:


plot(res$log2FoldChange, -log10(res$padj))

library(ggplot2)

# Assuming `res` is your DESeq2 results data frame
# Add a column for gene classification
res$threshold <- "NS"  # Not significant by default
res$threshold[res$padj < 0.05 & res$log2FoldChange > 0] <- "Upregulated in High-IL6 group"
res$threshold[res$padj < 0.05 & res$log2FoldChange < 0] <- "Downregulated in High-IL6 group"

# Convert to factor for consistent color mapping
res$threshold <- factor(res$threshold, levels = c("Upregulated in High-IL6 group", "Downregulated in High-IL6 group", "NS"))

# Volcano plot
volc <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = threshold)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("Upregulated in High-IL6 group" = "red", "Downregulated in High-IL6 group" = "blue", "NS" = "gray")) +
  xlim(-2, 2) +
  theme_minimal() +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Gene Regulation"
  ) +
  theme(legend.position = "right")

emf("Original_data/volcplot.emf", width = 6, height = 6)
# Your plot
volc
# Close device
dev.off()



resdf <- as.data.frame(res)


dfwbscapm$ex <- t(log10(count_eb["ENSG00000163568",]))
dfwbscapm$ex <- t(log10(count_eb["ENSG00000047410",]))

boxplot(dfwbscapm$ex ~ dfwbscapm$HiIL6)


upres <- resdf[resdf$stat > 0,]

head(upres[order(upres$pvalue),],20)




sigres <- resdf[!is.na(resdf$padj),]


library(org.Hs.eg.db)
sigres$entrez <-  mapIds(org.Hs.eg.db,
                         keys = row.names(sigres),
                         column = "ENTREZID",
                         keytype = "ENSEMBL",   # use ENSEMBL here
                         multiVals = "first")   # 


library(ReactomePA)


sigres  <- sigres[order(-sigres$stat),]

geneList <- sigres$stat

names(geneList) <- sigres$entrez

y <- gsePathway(geneList, 
                pvalueCutoff = 0.2,
                pAdjustMethod = "BH", 
                verbose = FALSE)
head(y, 30)


ydf <- as.data.frame(y)

RPR <- read.table("Original_data/reactomepathways.txt")





parents <- c("R-HSA-168256","R-HSA-1643685","R-HSA-109582")

length(unique(parents))

gen0 <- data.frame(ID=parents, gen=0)
gen1 <- data.frame(ID=subset(RPR, V1 %in% gen0$ID)$V2, gen=1)
gen2 <- data.frame(ID=subset(RPR, V1 %in% gen1$ID)$V2, gen=2)
gen3 <- data.frame(ID=subset(RPR, V1 %in% gen2$ID)$V2, gen=3)
gen4 <- data.frame(ID=subset(RPR, V1 %in% gen3$ID)$V2, gen=4)
gen5 <- data.frame(ID=subset(RPR, V1 %in% gen4$ID)$V2, gen=5)
gen6 <- data.frame(ID=subset(RPR, V1 %in% gen5$ID)$V2, gen=6)
gen7 <- data.frame(ID=subset(RPR, V1 %in% gen6$ID)$V2, gen=7)

allgens <- rbind(gen0,gen1,gen2,gen3,gen4,gen5,gen6,gen7)



yss <- ydf[ydf$ID %in% allgens$ID, ]


library(ggplot2)

# Ensure Description is a factor ordered by NES so bars are sorted
yss$Description <- factor(yss$Description, levels = yss$Description[order(yss$NES)])

# Plot
pathwayplot <- ggplot(yss, aes(x = NES, y = Description, fill = NES)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(x = "Normalized Enrichment Score (NES)", y = "Pathway") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10))



ggsave("Figure/pathwayplot.svg", plot = pathwayplot, width = 8, height = 10)

library(devEMF)

emf("Figure/pathwayplot.emf", width = 8, height = 10)
# Your plot
pathwayplot
# Close device
dev.off()
