rm(list = ls())

# Load required libraries
library(esvis)  # Load library for effect size visualization
library(pheatmap)  # Load library for creating heatmaps
library(tidyverse)  # Load library for data manipulation and visualization
library(dplyr)

#
Luminex_CAP_COVID_HV <- read.csv("Original_data/Luminex_CAP_COVID_HV_log.csv")

studydata_patients <- read.csv("Original_data/studydata_patients.csv")
studydata_HV <- read.csv("Original_data/studydata_HV.csv")
studydata <- rbind(studydata_patients, studydata_HV)

data <- merge(Luminex_CAP_COVID_HV, studydata, 
              by.x = "ID", by.y = "EB_id")

#compare medium between HV and CAP
data_medium <- data %>% filter(stimulation == "M")
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
  result <- hedg_g(data_medium, 
                   as.formula(paste(var, "~ group")), # column with group information
                   ref_group = "HV") # specify the reference group
  result <- mutate(result, Dependent_Variable = var)
  results_list[[var]] <- result
}

# Combine the results into a single DataFrame
combined_results <- bind_rows(results_list)

# Print the combined DataFrame
print(combined_results)

# Prepare data for the heatmap
data_heatmap_hedg_g <- dplyr::select(combined_results, Dependent_Variable, hedg_g, group_foc)  # Select relevant columns
data_heatmap_hedg_g <- data_heatmap_hedg_g %>%
  pivot_wider(names_from = group_foc, values_from = hedg_g)  # Reshape data for heatmap
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
output_pdf_file <- "Figure/variables_pro_CAP_HV.pdf"
pdf(file = output_pdf_file, width = 2, height = 3)


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


#####
data_sti <- data %>% filter(!stimulation == "M")
# List of domains
variables_chemo <- c("minis_M_CCL2","minis_M_CCL3","minis_M_CCL4")
variables_ati <- c("minis_M_IL_8_1","minis_M_TNF","minis_M_IL_6","minis_M_IL_1beta")
variables_pro <- c("minis_M_IL_1RA","minis_M_IL_10")

# List to store the results
results_list <- list()
# Perform hedg_g for each dependent variable
####if we want to get the Hedges'g for different domains, we need to change the variables_domains##############################################################
for (var in variables_chemo) {
  # Conduct hedg_g analysis
  result <- hedg_g(data_sti, 
                   as.formula(paste(var, "~ group")), # column with group information
                   ref_group = "HV") # specify the reference group
  result <- mutate(result, Dependent_Variable = var)
  results_list[[var]] <- result
}

# Combine the results into a single DataFrame
combined_results <- bind_rows(results_list)

# Print the combined DataFrame
print(combined_results)

# Prepare data for the heatmap
data_heatmap_hedg_g <- dplyr::select(combined_results, Dependent_Variable, hedg_g, group_foc)  # Select relevant columns
data_heatmap_hedg_g <- data_heatmap_hedg_g %>%
  pivot_wider(names_from = group_foc, values_from = hedg_g)  # Reshape data for heatmap
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
output_pdf_file <- "Figure/variables_chemo_minis_CAP_HV.pdf"
pdf(file = output_pdf_file, width = 2, height = 3)


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