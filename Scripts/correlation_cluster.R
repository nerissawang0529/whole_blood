### make matrix of piD versus "markers" (cytokine_stimulation combinations n=36)
rm(list = ls())
pacman::p_load(pacman, dplyr, readr, edgeR, AnnotationDbi, org.Hs.eg.db)

studydata_patients <- read.csv("Original_data/studydata_patients.csv")
#studydata_HV <- read.csv("Original_data/studydata_HV.csv")
#studydata_combined <- rbind(studydata_patients, studydata_HV)

Luminex_CAP_COVID_HV <- read.csv("Original_data/Luminex_CAP_COVID_HV.csv")
merged_data <- merge(Luminex_CAP_COVID_HV, studydata_patients, 
                     by.x = "ID", by.y = "EB_id")

mdata <- merged_data[,1:12]
mdata <- mdata %>% 
  filter(stimulation != "M")

mdata$Day <- NULL

library(reshape2)
ldata <- reshape2::melt(mdata, id=c("ID","stimulation") )
ldata$level <- log10(ldata$value)

ldata$value <- NULL


raw_w <- reshape2::dcast(ldata, ID ~ stimulation + variable, value.var="level" )

row.names(raw_w) <- raw_w$ID
raw_w$ID <- NULL

corm <- cor(raw_w, method = "spearman",
            use = "pairwise.complete.obs")

## After you compute `corm` ##

# helper: map stim prefixes
stim_map <- c(
  "LPS" = "LPS",
  "KP"  = "Kp",   # Klebsiella pneumoniae
  "PP"  = "Spn"   # Streptococcus pneumoniae (pneumococcus)
)

# helper: map cytokine tokens to nicer labels
marker_map <- c(
  "IL_1beta" = "IL-1\u03B2",  # IL-1β
  "IL_1RA"   = "IL-1RA",
  "IL_8_1"   = "IL-8",
  "IL_6"     = "IL-6",
  "IL_10"    = "IL-10",
  "TNF"      = "TNF",
  "CCL2"     = "CCL2",
  "CCL3"     = "CCL3",
  "CCL4"     = "CCL4"
)

pretty_names <- function(nms){
  # split "KP_IL_6" -> stim="KP", mark="IL_6"
  stim <- sub("^([A-Za-z]+)_.*$", "\\1", nms)
  mark <- sub("^[A-Za-z]+_(.*)$", "\\1", nms)
  
  stim_lab <- ifelse(stim %in% names(stim_map), stim_map[stim], stim)
  mark_lab <- ifelse(mark %in% names(marker_map), marker_map[mark], mark)
  
  paste0(stim_lab, "–", mark_lab)   # e.g., "Kp–IL-6"
}

new_names <- pretty_names(colnames(corm))
colnames(corm) <- rownames(corm) <- new_names




library(corrplot)
corrplot(corm)
corrplot(corm, method = "color", type = "upper", 
         col = colorRampPalette(c("blue", "white", "red"))(200),
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
        color = colorRampPalette(c("blue", "white", "red"))(100),
        breaks = breaks,
        main = "Spearman Correlation Heatmap")


nclust <- cutree(hco, k=3)
table(nclust)

# Convert to a data frame for easy viewing
cluster_assignments <- data.frame(Marker = names(nclust), Cluster = nclust)

# View which markers belong to which cluster
cluster_assignments <- cluster_assignments %>%
  arrange(Cluster)
print(cluster_assignments)



# distance + clustering (you already have these)
D   <- as.dist(1 - corm)
hco <- hclust(D, method = "complete")

library(corrplot)
corrplot(
  corm,
  method = "color",
  type   = "upper",            # <- only half
  order  = "hclust",
  hclust.method = "complete",
  col = colorRampPalette(c("blue","white","red"))(200),
  tl.col = "black", tl.srt = 60,
  diag = FALSE                 # hide diagonal
)


## ---- trait parser on your column names like "KP_IL_6" ----
trait_df <- function(mat) {
  tibble::tibble(
    trait  = colnames(mat),
    stim   = sub("^([A-Za-z]+)_.*$", "\\1", colnames(mat)),
    marker = sub("^[A-Za-z]+_(.*)$", "\\1", colnames(mat))
  )
}

## ---- pretty labeler (uses your maps) ----
pretty_trait <- function(trait_vec){
  stim <- sub("^([A-Za-z]+)_.*$", "\\1", trait_vec)
  mark <- sub("^[A-Za-z]+_(.*)$", "\\1", trait_vec)
  stim_lab <- ifelse(stim %in% names(stim_map), stim_map[stim], stim)
  mark_lab <- ifelse(mark %in% names(marker_map), marker_map[mark], mark)
  paste0(stim_lab, "–", mark_lab)
}

## ---- generic: subset columns, compute Spearman, plot heatmap ----
make_heatmap <- function(mat, cols, title = "Spearman correlation heatmap") {
  stopifnot(length(cols) >= 2)
  submat <- mat[, cols, drop = FALSE]
  R <- cor(submat, use = "pairwise.complete.obs", method = "spearman")
  colnames(R) <- rownames(R) <- pretty_trait(colnames(R))
  D <- as.dist(1 - R)
  pheatmap::pheatmap(
    R,
    clustering_distance_rows = D,
    clustering_distance_cols = D,
    color  = colorRampPalette(c("blue","white","red"))(100),
    breaks = seq(-1, 1, length.out = 101),
    main   = title
  )
  invisible(R)
}

## define groups based on your raw marker tokens
innate_set    <- c("IL_1beta","IL_6","TNF","IL_10","IL_8_1","IL_1RA")
chemokine_set <- c("CCL2","CCL3","CCL4")

td <- trait_df(raw_w)

## Innate only (all stimuli kept)
cols_innate <- td$trait[td$marker %in% innate_set]
R_innate <- make_heatmap(raw_w, cols_innate,
                         title = "Innate cytokines (all stimuli)")

## Chemokines only (all stimuli kept)
cols_chem <- td$trait[td$marker %in% chemokine_set]
R_chem <- make_heatmap(raw_w, cols_chem,
                       title = "Chemokines (all stimuli)")

## choose the stimuli you have
stims <- c("LPS","KP","PP")  # adjust if needed

for (s in stims) {
  cols_s <- trait_df(raw_w) |>
    dplyr::filter(stim == s) |>
    dplyr::pull(trait)
  if (length(cols_s) >= 2) {
    make_heatmap(raw_w, cols_s,
                 title = paste0("Within-stimulus module: ", s))
  }
}


## pick a marker token as it appears in your column names
marker_pick <- "IL_6"   # try "TNF", "IL_1beta", "CCL3", etc.

cols_m <- trait_df(raw_w) |>
  dplyr::filter(marker == marker_pick) |>
  dplyr::pull(trait)

if (length(cols_m) >= 2) {
  ttl <- paste0("Across-stimulus module for ", marker_map[[marker_pick]] %||% marker_pick)
  make_heatmap(raw_w, cols_m, title = ttl)
}

