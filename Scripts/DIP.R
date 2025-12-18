# =================== 0. Install and set up Python virtualenv ===================
library(reticulate)

# Create virtualenv for DIP if not already created
if (!virtualenv_exists("dip-env")) {
  virtualenv_create("dip-env")
}

# Activate the virtualenv
use_virtualenv("dip-env", required = TRUE)

# Install required Python modules into the virtualenv
py_install(c("scikit-learn", "numpy", "pandas"), pip = TRUE)

# Confirm scikit-learn is available
stopifnot(py_module_available("sklearn"))

# =================== 1. Install DIP package from GitHub (using pak) ===================
# If not yet installed, run once:
install.packages("pak")  # Safe alternative to devtools

# Use pak to install the DIP package from GitHub
pak::pkg_install("DysregulatedImmuneProfile/DIP")

# Load DIP
library(DIP)

# =================== 2. Load required packages ===================
# Use pacman to load/install R packages
install.packages("pacman")
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr,
               Hmisc, tableone, kmed, survival, survminer, skimr, grDevices, ggrepel)

# =================== 3. Load and filter data ===================
# Read plasma biomarker data
direct_plasma_patients <- read.csv("Original_data/first_syndecan_merged_markers_raw.csv")

# Filter patients with ID containing 'EB'
direct_plasma_EBpatients <- direct_plasma_patients %>%
  filter(str_detect(New_ID, "EB")) %>%
  mutate(ID = str_extract(New_ID, "\\d+$"))

# Check overlap with luminex data (luminex_long must be loaded earlier)
table(direct_plasma_EBpatients$ID %in% luminex_long$ID)
table(unique(luminex_long$ID) %in% direct_plasma_EBpatients$ID)

# Filter to patients with matched stimulation data
filtered_EBpatients <- direct_plasma_EBpatients %>%
  filter(ID %in% unique(luminex_long$ID))

# Select 3 markers needed for DIP score
selected_EBpatients <- filtered_EBpatients %>%
  select(ID, TREM_1, IL_6, Procalcitonin)

# =================== 4. Export data to CSV ===================
destination_folder <- "Original_data"
export_file_name <- "DIP_file.csv"
write.csv(selected_EBpatients, file.path(destination_folder, export_file_name), row.names = FALSE)
cat("File exported to:", file.path(destination_folder, export_file_name), "\n")

# =================== 5. Run DIP Shiny App ===================
DIP_app()  # This should now run without Python error

# =================== 6. (Optional) Directly compute cDIP scores ===================
test_data <- data.frame(
  TREM_1 = c(100, 200, 300),
  IL_6 = c(500, 1000, 1500),
  Procalcitonin = c(50, 100, 150)
)

cDIP(test_data)  # This should return DIP scores
