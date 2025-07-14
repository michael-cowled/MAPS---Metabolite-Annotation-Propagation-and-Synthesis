##Filename extractor for GNPS metadata creation

extract_filenames_from_subfolders <- function(directory_path) {
  # Check if the directory exists
  if (!dir.exists(directory_path)) {
    stop("Error: The specified directory does not exist.")
  }
  
  # Get all subdirectories within the main directory
  subfolders <- list.dirs(directory_path, full.names = TRUE, recursive = FALSE)
  
  # Initialize an empty list to store data frames for each subfolder
  all_files_list <- list()
  
  # Loop through each subfolder
  for (subfolder_path in subfolders) {
    # Extract the subfolder name
    subfolder_name <- basename(subfolder_path)
    
    # Get all files within the current subfolder
    files_in_subfolder <- list.files(subfolder_path, full.names = FALSE)
    
    # If there are files in the subfolder, create a data frame for them
    if (length(files_in_subfolder) > 0) {
      temp_df <- data.frame(
        filename = files_in_subfolder,
        subfolder = rep(subfolder_name, length(files_in_subfolder)),
        stringsAsFactors = FALSE
      )
      all_files_list[[subfolder_name]] <- temp_df
    }
  }
  
  # Combine all individual data frames into one
  if (length(all_files_list) > 0) {
    result_df <- do.call(rbind, all_files_list)
    row.names(result_df) <- NULL # Reset row names
    return(result_df)
  } else {
    message("No subfolders containing files were found in the specified directory.")
    return(data.frame(filename = character(0), subfolder = character(0)))
  }
}

main_dir <- "Y:/MA_BPA_Microbiome/LCMS data/HGMD_0129_ReDU/ReDU"

 df <- extract_filenames_from_subfolders(main_dir)
