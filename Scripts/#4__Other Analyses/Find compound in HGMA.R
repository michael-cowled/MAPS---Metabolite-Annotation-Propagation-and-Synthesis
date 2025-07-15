# Load necessary libraries
library(dplyr)
library(readr)
library(purrr) # Already there, but useful for more advanced functional programming
library(tidyr) # For pivot_longer, pivot_wider etc.
library(tidyverse) # Loads a collection of packages including dplyr, readr, tidyr, purrr, stringr, ggplot2, forcats

# --- Load Data ---
# It's good practice to wrap file paths in tryCatch to handle potential errors
# if the files are not found.

tryCatch({
  # Retrieve full sample metadata
  metadata <- read.csv("C:/Users/mcowled/The University of Melbourne/MA Human Gut Metabolome - Documents/Data Management and Analysis/HGM_data_merge_metadata/full_metadata_prefixed_joined_simple_join.csv")
  
  # Retrieve all LCMS data
  all_merged_df <- read.csv("C:/Users/mcowled/The University of Melbourne/MA Human Gut Metabolome - Documents/Data Management and Analysis/HGM_data_merge_metadata/ALL_MERGED_SAMPLE_ANNOTATION.csv") %>%
    filter(Best.Annotation.Confidence.Level <= 3)
  
  # Retrieve all hits from immune screen
  hits <- read.csv("C:/Users/mcowled/The University of Melbourne/MA Human Gut Metabolome - Documents/Data Management and Analysis/Standards and Compounds Australia/PMC_266_CompAust_Hits.csv")
  
  CA.lib <- read.csv("C:/Users/mcowled/The University of Melbourne/MA Human Gut Metabolome - Documents/Data Management and Analysis/Standards and Compounds Australia/CompAus_Overlapping_CIDs_with_fecal_fractions.csv")
  
  Best.Annotations <- read.csv("Y:/MA_BPA_Microbiome/Total-List-Of-Annotations.csv")
  
}, error = function(e) {
  stop("Error loading one or more CSV files. Please check the file paths. Original error: ", e$message)
})

# ------------------------------------------------------------------------------#
# Prepare hits and all_merged_df for matching
# This block is critical for setting up the CID columns correctly before the loop.

# 1. Prepare CA.lib and hits: Add CID to hits based on Compound_Name from CA.lib
CA.lib <- CA.lib %>%
  rename(Compound_Name = Product.Name)

# Perform a left join to add 'cid' from CA.lib to 'hits'.
# Ensure 'cid' column from CA.lib is named 'cid' in hits after this join.
hits <- hits %>%
  left_join(CA.lib %>% select(Compound_Name, cid), by = "Compound_Name")

# 2. Prepare all_merged_df: Add CID from Best.Annotations
# This ensures 'all_merged_df' has a 'CID' column to be used for fallback matching.
all_merged_df <- all_merged_df %>%
  left_join(Best.Annotations %>% select(Best.Annotation, CID),
            by = "Best.Annotation")

# Rename the 'CID' column in all_merged_df for consistency if needed.
# Your original script renames 'cid' to 'CID' but 'all_merged_df' gets 'CID' directly from Best.Annotations.
# This line might be redundant or could cause issues if 'cid' was not pre-existing in all_merged_df.
# If Best.Annotations.CID is already named 'CID', this line is not needed.
# If it was named something else (e.g., Best.Annotations.CID_Value), you'd rename it here.
# Assuming Best.Annotations$CID correctly provides 'CID', no rename needed here for 'cid' to 'CID'.
# all_merged_df <- all_merged_df %>% rename(CID = cid) # <--- Commenting out unless 'cid' truly exists and needs renaming to 'CID'

# Now, ensure that 'hits' has a 'cid' column (from CA.lib join) and 'all_merged_df' has a 'CID' column (from Best.Annotations join)

# --- Add a unique identifier to 'hits' for easier joining in the loop ---
# This is crucial for handling multiple matches or to track original hit rows.
hits <- hits %>%
  mutate(HitID = row_number()) # Add a unique ID to each row in hits

# -------------------------------------------------------------------------------#

# --- Data Preparation ---

# Replace NA in 'area' with 0 for consistent filtering
all_merged_df <- all_merged_df %>%
  mutate(area = ifelse(is.na(area), 0, area))

# Join the LCMS data with the sample metadata
# Note: This can be a large object. Ensure you have enough memory.
all_merged_df_with_metadata <- all_merged_df %>%
  left_join(metadata, by = c("samples" = "A.HGMA.ID"))


# --- Loop Through Each Compound and Export Data ---

# Get a unique list of HitIDs from the 'hits' dataframe
# We'll iterate by HitID to correctly apply fallback logic per original hit.
# This ensures each entry in 'hits' is processed individually, enabling the fallback.
hit_ids_to_process <- unique(hits$HitID)

# Create a directory to store the output files if it doesn't exist
output_dir <- "compound_queries"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Loop over each HitID from the 'hits' dataframe
for (current_hit_id in hit_ids_to_process) {
  
  # Extract the specific hit row for the current iteration
  current_hit <- hits %>%
    filter(HitID == current_hit_id)
  
  # Get the Compound_Name and CID from the current hit
  target_compound_name <- current_hit$Compound_Name
  target_cid_from_hit <- current_hit$cid # This is the CID from hits (from CA.lib)
  
  cat(paste("Processing Hit ID:", current_hit_id, "- Compound:", target_compound_name, "(Fallback CID:", target_cid_from_hit, ")\n"))
  
  # --- Primary Match: Filter by Compound_Name ---
  target_rows_by_compound <- all_merged_df_with_metadata %>%
    filter(!is.na(area), area > 0) %>%
    filter(Best.Annotation == target_compound_name)
  
  # Check if primary match found data
  if (nrow(target_rows_by_compound) > 0) {
    target_rows <- target_rows_by_compound
    match_type <- "by_Compound_Name"
    cat(paste("  -> Found", nrow(target_rows), "rows by Compound_Name.\n"))
  } else {
    # --- Fallback Match: If Compound_Name fails, try matching by CID (from hits) ---
    if (!is.na(target_cid_from_hit)) { # Only try fallback if a CID exists in hits
      target_rows_by_cid <- all_merged_df_with_metadata %>%
        filter(!is.na(area), area > 0) %>%
        filter(CID == target_cid_from_hit) # Match with 'CID' in all_merged_df_with_metadata
      
      if (nrow(target_rows_by_cid) > 0) {
        target_rows <- target_rows_by_cid
        match_type <- "by_CID_fallback"
        cat(paste("  -> No match by Compound_Name. Found", nrow(target_rows), "rows by CID (fallback).\n"))
      } else {
        target_rows <- data.frame() # No rows found
        match_type <- "none"
        cat(paste("  -> No matching data found by Compound_Name or CID for", target_compound_name, "\n"))
      }
    } else {
      target_rows <- data.frame() # No rows found and no CID for fallback
      match_type <- "none"
      cat(paste("  -> No matching data found by Compound_Name and no CID available for fallback for", target_compound_name, "\n"))
    }
  }
  
  # --- Export Data (if any rows found) ---
  if (nrow(target_rows) > 0) {
    # Sanitize the compound name to create a valid filename
    # Combine compound name and HitID for a unique filename
    safe_filename <- gsub("[^a-zA-Z0-9_]", "_", paste0(target_compound_name, "_HitID_", current_hit_id, "_", match_type))
    
    # Define the output file path
    output_path <- file.path(output_dir, paste0("compound_query_", safe_filename, ".csv"))
    
    # Write the filtered data to a CSV file
    write.csv(target_rows, output_path, row.names = FALSE)
    
    cat(paste("  -> Data saved to:", output_path, "\n\n"))
    
  } else {
    cat(paste("  -> No data exported for Hit ID", current_hit_id, "\n\n"))
  }
}

cat("Processing complete.\n")