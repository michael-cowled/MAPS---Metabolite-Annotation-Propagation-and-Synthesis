library(dplyr)
library(readr)

#This chunk is to plot the count of unique annotations in each dataset eg. "HGMD_0110" compared to other datasets aaccording to the superclass of the annotations.
# Need to run this chunk after the Comparing_datasets_MC.R
# set the working directory to the working folder for dataset comparison, CHANGE THE FILDER NAME BEFORE RUNNING!
setwd("C:\Users\mcowled\Documents\HGMD_0075_vs_0108")

# Load filtered_datasets
filtered_datasets <- read_csv("filtered_datasets.csv")

# Define base path for original annotation files
base_path <- "Y:/MA_BPA_Microbiome/Dataset-Annotations"

# Add empty Superclass column
filtered_datasets$Superclass <- NA

# Function to choose the correct superclass based on annotation type
get_superclass <- function(type, gnps, ms2query, canopus) {
  if (is.na(type) || type %in% c("csi", "canopus", "Propagated Analogue")) {
    return(canopus)
  } else if (type %in% c("ms2query", "ms2query analogue")) {
    return(ms2query)
  } else if (type == "gnps") {
    return(gnps)
  } else {
    return(NA)
  }
}

# Loop through each unique dataset.id
for (dataset in unique(filtered_datasets$dataset.id)) {
  # Path to the corresponding original dataset
  dataset_path <- file.path(base_path, paste0(dataset, ".csv"))
  
  # Load the original dataset
  original_data <- read_csv(dataset_path, show_col_types = FALSE)
  
  # Rename the first column to 'Feature_ID' if it's unnamed
  if (names(original_data)[1] == "...1") {
    names(original_data)[1] <- "Feature_ID"
  }
  
  # Filter rows in filtered_datasets for the current dataset
  idx <- filtered_datasets$dataset.id == dataset
  
  # Join filtered_datasets with original_data by Best.Annotation.Smiles
  merged <- filtered_datasets[idx, ] %>%
    left_join(original_data, by = "Best.Annotation.Smiles")
  
  # Get superclass values using the defined function
  merged$Superclass <- mapply(
    get_superclass,
    merged$Best.Annotation.Type,
    merged$gnps.NPC.superclass,
    merged$ms2query.NPC.superclass,
    merged$canopus.NPC.superclass
  )