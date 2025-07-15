##File List Extractor##

# Specify the path to the Data Management Plan (unless moved to Mediaflux then should be consistent)
excel_file <- "C:/Users/mcowled/The University of Melbourne/Gut microbiome group - Documents/Data Management and Analysis/HGM - Data Management System.xlsx"

## Loading libraries

# Function to check and install missing packages
check_and_install <- function(packages) {
  # Identify packages that are not installed
  to_install <- packages[!(packages %in% installed.packages()[, "Package"])]
  # Install missing packages
  if (length(to_install) > 0) {
    install.packages(to_install)
  }  # Load all required packages
  invisible(lapply(packages, library, character.only = TRUE))
}
required_packages <- c("dplyr", "tidyr", "stringr", "readr", "reshape2", "ggplot2", "svglite", "readxl", "openxlsx") # Check, install, and load required packages
check_and_install(required_packages) # Now your packages are installed and loaded

## Updating Metadata                                                              ----> Remove from published version

sheet_names <- excel_sheets(excel_file) # Get the sheet names
for (sheet in sheet_names) {
  data <- read_excel(excel_file, sheet = sheet)
  write.csv(data, file = paste0(sheet, ".csv"), row.names = FALSE)
}

# Dataset ID: From Data Management Plan:
dataset <- read.csv("HGM/D - Dataset.csv")
analyses <- read.csv("HGM/A - Analysis.csv") # Read analyses file *before* the loop

# Get unique HGMD.ID values
hgmd_ids <- unique(dataset$HGMD.ID)

# Iterate over unique HGMD.ID values
for (dataset.id in hgmd_ids) {
  # Filter dataset for the current HGMD.ID
  current_dataset <- dataset %>%
    filter(HGMD.ID == dataset.id)
  
  # Check if any matching HGMD.ID was found. If not, skip to the next ID.
  if (nrow(current_dataset) == 0) {
    print(paste("No data found for HGMD.ID:", dataset.id))
    next # Skip to next iteration
  }
  
  output_directory <- current_dataset$Processed.Data.Folder[1]
  mzml.folder <- paste0(current_dataset$Processed.Data.Folder, "\\mzml")
  data.folder <- paste0(current_dataset$Processed.Data.Folder, "\\data")
  print(folder)
  
  # Define the function to process the folder and files (same as before)
  process_files <- function(mzml.folder, data.folder) {
    if (dir.exists(mzml.folder)) {
      files <- list.files(mzml.folder)
      full_paths <- file.path(mzml.folder, files)
      full_paths_windows <- gsub("/", "\\\\", full_paths)
      file_info <- data.frame(
        mzml.folder = mzml.folder,
        filename = files,
        full_path = full_paths_windows
      )
      return(file_info)
    } else if (dir.exists(data.folder)) {
      files <- list.files(data.folder)
      full_paths <- file.path(data.folder, files)
      full_paths_windows <- gsub("/", "\\\\", full_paths)
      file_info <- data.frame(
        data.folder = data.folder,
        filename = files,
        full_path = full_paths_windows
      )
      return(file_info)
    } else {
      print(paste("Warning: Folder does not exist:", mzml.folder)) # Print a warning, don't stop the loop
      return(data.frame()) # Return an empty dataframe to avoid errors later
    }
  }
  
  file_info <- process_files(mzml.folder, data.folder)
  
  # Filter analyses using *string contains* logic for dataset.ID
  analyses2 <- analyses %>%
    filter(grepl(dataset.id, dataset.ID)) # Changed to grepl for "contains" logic
  
  if (nrow(analyses2) > 0 && nrow(file_info) > 0) {
    for (i in 1:nrow(analyses2)) {
      
      filtered_df <- file_info[grepl(analyses2[i, 1], file_info$filename), ]
      
      if (nrow(filtered_df) > 0) {
        print(filtered_df[1, 3])
        file_path <- filtered_df[1, 3]
        
        # Handle multiple matching rows in 'analyses'
        matching_row_indices <- which(analyses$HGMA.ID == analyses2[i, 1]) #Find all matching rows
        
        for (j in matching_row_indices) { # Loop through all matched rows.
          analyses$File.Location[j] <- gsub("\\\\", "/", file_path) #Update all the matching rows
        }
        
        
      } else {
        print(paste("Warning: No matching files found for HGMA.ID:", analyses2[i, 1], "in folder:", folder))
      }
    }
  } else {
    print(paste("Warning: No matching analyses or files found for dataset.ID (containing):", dataset.id))
  }
}


##Phase 3. Update to Excel
wb <- loadWorkbook(excel_file)
writeData(wb, sheet = "A - Analysis", analyses)
saveWorkbook(wb, excel_file, overwrite = TRUE)