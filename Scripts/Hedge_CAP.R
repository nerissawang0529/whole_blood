rm(list = ls())

# Load required libraries
library(esvis)  # Load library for effect size visualization
library(pheatmap)  # Load library for creating heatmaps
library(tidyverse)  # Load library for data manipulation and visualization
library(dplyr)

#
Luminex_CAP_COVID_HV <- read.csv("Original_data/Luminex_CAP_COVID_HV_log.csv")
studydata_patients <- read.csv("Original_data/studydata_patients.csv")
data <- merge(Luminex_CAP_COVID_HV, studydata_patients, 
                     by.x = "ID", by.y = "EB_id")


# List of domains
variables_chemo <- c("CCL2","CCL3","CCL4")
variables_ati <- c("IL_8_1","TNF","IL_6","IL_1beta")
variables_pro <- c("IL_1RA","IL_10")

# List to store the results
results_list <- list()
# Perform hedg_g for each dependent variable
####if we want to get the Hedges'g for different domains, we need to change the variables_domains##############################################################
for (var in variables_pro) {
  # Conduct hedg_g analysis
  result <- hedg_g(data, 
                   as.formula(paste(var, "~ stimulation")), # column with group information
                   ref_group = "M") # specify the reference group
  result <- mutate(result, Dependent_Variable = var)
  results_list[[var]] <- result
}

# Combine the results into a single DataFrame
combined_results <- bind_rows(results_list)

# Print the combined DataFrame
print(combined_results)

# Prepare data for the heatmap
data_heatmap_hedg_g <- dplyr::select(combined_results, Dependent_Variable, hedg_g, stimulation_foc)  # Select relevant columns
data_heatmap_hedg_g <- data_heatmap_hedg_g %>%
  pivot_wider(names_from = stimulation_foc, values_from = hedg_g)  # Reshape data for heatmap
data_heatmap_hedg_g <- column_to_rownames(data_heatmap_hedg_g, var = "Dependent_Variable")  # Set variable as row names


# Print the heatmap data
print(data_heatmap_hedg_g)


# Define breaks and colors for the heatmap
breaks <- c(-Inf, -1.5, -0.8, -0.5, -0.2, 0.2, 0.5, 0.8, 1.5, Inf)


#for the figure
##Question:Is the color I set correct?
colors <- c("#f70d1a","#ff533f","#ff9e88","#ffbfae", "white", "#9F9FFF","#4444FF","#0000E7","#00008B")





####if we want to get the Hedges'g for different domains, we need to change the variables_domains##############################################################
#output the figure
output_pdf_file <- "Figure/variables_pro_CAP.pdf"
pdf(file = output_pdf_file, width = 3, height = 2)


# Create the heatmap using pheatmap with custom breaks and colors
pheatmap(data_heatmap_hedg_g,
         scale = "none",
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize = 12,
         border_color = "black",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         width = 2, height = 5, dpi = 600,
         color = colors,
         breaks = breaks,
         legend = FALSE
)

dev.off()

