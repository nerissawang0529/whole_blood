### summary TNF

ldata_sel <- ldata[ldata$stimulation %in% c("PP","LPS","KP") &
                     ldata$variable %in% c("TNF"),]

ldata_wide_ss <- reshape2::dcast(ldata_sel, ID ~ stimulation, value.var = "level")

row.names(ldata_wide_ss) <- ldata_wide_ss$ID
ldata_wide_ss$ID <- NULL


pca_ss <- prcomp(ldata_wide_ss)

library(ggfortify)

autoplot(pca_ss, loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3 )


pca_df <- data.frame(pca_ss$x[, 1:2])
pca_df$ID <- as.integer(rownames(pca_df))
pca_df$group <- "CAP"
pca_df$cluster <- "CAP"

# project healthies using predict

ldata_HV_sel <- ldata_HV[ldata_HV$stimulation %in% c("PP","LPS","KP") &
                           ldata_HV$variable %in% c("TNF"),]

ldata_HV_wide_ss <- reshape2::dcast(ldata_HV_sel, ID ~ stimulation, value.var = "level")

row.names(ldata_HV_wide_ss) <- ldata_HV_wide_ss$ID
ldata_HV_wide_ss$ID <- NULL


pca_HV <- predict(pca_ss, newdata = ldata_HV_wide_ss)



pca_HV_df <- data.frame(pca_HV[, 1:2])
pca_HV_df$ID <- as.integer(rownames(pca_HV_df))
pca_HV_df$group <- "Healthy group"
pca_HV_df$cluster <- "HV"

pca_all <- rbind(pca_df, pca_HV_df)

# ==== 8. Loadings (arrow) ====
library(tidyverse)
loadings_df <- as.data.frame(pca_ss$rotation[, 1:2])
loadings_df$variable <- rownames(loadings_df)
loadings_df$stim <- sapply(strsplit(loadings_df$variable, "_"), function(parts) tail(parts, 1))
loadings_df$marker <- str_remove( rownames(loadings_df), "_[^_]+$")
loadings_df$label <- loadings_df$marker

# ==== 9. PCA plot ====
arrow_scale <- 6
label_offset <- 1.1
label_size <- 5
marker_colors <- c( "TNF" = "#FF7F00")

library(ggrepel)

ggplot(pca_all, aes(PC1, PC2, color = cluster, shape = group)) +
  geom_point(size = 2, alpha = 0.85) +
  stat_ellipse(aes(group = cluster), linewidth = 1.2) +
  geom_segment(data = loadings_df,
               aes(x = 0, y = 0,
                   xend = PC1 * arrow_scale,
                   yend = PC2 * arrow_scale,
                   color = stim),
               arrow = arrow(length = unit(0.35, "cm")),
               linewidth = 1.2,
               inherit.aes = FALSE) +
  geom_label_repel(data = loadings_df,
                   aes(x = PC1 * arrow_scale * label_offset,
                       y = PC2 * arrow_scale * label_offset,
                       label = label, fill = marker),
                   color = "black",
                   size = label_size,
                   label.padding = 0.3,
                   label.size = 0.4,
                   segment.color = "grey30",
                   max.overlaps = Inf,
                   box.padding = 0.5,
                   inherit.aes = FALSE) 


boxplot(pca_all$PC1 ~ pca_all$group)
