rm(list = ls())

# Read MS1 gene expression matrix from Reyes et al.
ms1_raw <- read.table(
  "Original_data/MS1_septiStates_150_Reyes_2020.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# Extract MS1-related gene list (gene names only)
ms1_genes <- unique(ms1_raw$gennames)

library(AnnotationDbi)
library(reactome.db)

# Reactome ID for neutrophil degranulation
reactome_id <- "R-HSA-6798695"

neutrophil_reactome <- AnnotationDbi::select(
  reactome.db,
  keys     = reactome_id,
  keytype  = "PATHID",
  columns  = c("ENTREZID")
)

# Clean up
neutrophil_entrez <- unique(neutrophil_reactome$ENTREZID)
neutrophil_entrez <- neutrophil_entrez[!is.na(neutrophil_entrez)]


library(org.Hs.eg.db)

neutrophil_symbols <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys     = neutrophil_entrez,
  keytype  = "ENTREZID",
  columns  = c("SYMBOL")
)

neutrophil_genes <- unique(neutrophil_symbols$SYMBOL)
neutrophil_genes <- neutrophil_genes[!is.na(neutrophil_genes)]


# Compute overlap between MS1 and neutrophil gene sets
overlap_genes <- intersect(ms1_genes, neutrophil_genes)

# Summary
length(ms1_genes)
length(neutrophil_genes)
length(overlap_genes)



install.packages("VennDiagram")

library(VennDiagram)
library(grid)
gene_sets <- list(
  "MS1 signature" = ms1_genes,
  "Neutrophil degranulation" = neutrophil_genes
)
venn.plot <- venn.diagram(
  x = gene_sets,
  filename = NULL,          # draw to R device
  fill = c("#E69F00", "#56B4E9"),
  alpha = 0.6,
  cex = 1.4,                # number size
  cat.cex = 1.2,            # label size
  cat.pos = c(-20, 20),
  margin = 0.1
)

grid.newpage()
grid.draw(venn.plot)




# Read MS1 gene expression matrix from Reyes et al.
il1r2_raw <- read.table(
  "Original_data/CIBERSORTx_Job16_ref_matrix_equal_inferred_phenoclasses.CIBERSORTx_Job16_ref_matrix_equal_inferred_refsample.bm.K999.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)
# Extract IL1R2+ associated genes
il1r2_genes <- il1r2_raw$NAME[
  il1r2_raw$IL1R2._immature_neutrophils > 1
]

length(il1r2_genes)
length(neutrophil_genes)
overlap_il1r2_neutrophil <- intersect(il1r2_genes, neutrophil_genes)

length(overlap_il1r2_neutrophil)

gene_sets <- list(
  "IL1R2+ signature" = il1r2_genes,
  "Neutrophil degranulation" = neutrophil_genes
)
venn.plot <- venn.diagram(
  x = gene_sets,
  filename = NULL,          # draw to R device
  fill = c("#E69", "#56B4E9"),
  alpha = 0.6,
  cex = 1.4,                # number size
  cat.cex = 1.2,            # label size
  cat.pos = c(-20, 20),
  margin = 0.1
)

grid.newpage()
grid.draw(venn.plot)


# Overlap between MS1 and IL1R2+
overlap_ms1_il1r2 <- intersect(ms1_genes, il1r2_genes)
cat("MS1 genes:", length(ms1_genes), "\n")
cat("IL1R2+ genes:", length(il1r2_genes), "\n")
cat("Overlap:", length(overlap_ms1_il1r2), "\n")
il1r2_genes <- unique(as.character(il1r2_genes))
neutrophil_genes <- unique(as.character(neutrophil_genes))


gene_sets <- list(
  "IL1R2+ signature" = il1r2_genes,
  "MS1 signature" = ms1_genes
)

library(VennDiagram)
library(grid)

venn.plot <- venn.diagram(
  x = gene_sets,
  filename = NULL,          # draw to R device
  fill = c("#E69", "#E69F00"),
  alpha = 0.6,
  cex = 1.4,                # number size
  cat.cex = 1.2,            # label size
  cat.pos = c(-20, 20),
  cat.dist = c(0.05, 0.05), # label distance
  margin = 0.1
)

grid.newpage()
grid.draw(venn.plot)
