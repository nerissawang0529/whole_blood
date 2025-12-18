cap_healthy <- read.csv("Original_data/COCONUT_all_gene_expression_data.csv")
EB_all <- colnames(cap_healthy)[
  grepl("^X\\d+", colnames(cap_healthy)) &
    !grepl("OPT|MARS", colnames(cap_healthy))
]

length(EB_all)
head(EB_all, 20)
