# Load necessary libraries
library(dplyr)
library(readr)
library(purrr) # Already there, but useful for more advanced functional programming
library(tidyr) # For pivot_longer, pivot_wider etc.
library(tidyverse) # Loads a collection of packages including dplyr, readr, tidyr, purrr, stringr, ggplot2, forcats

# --- Load Data ---
# Load necessary library

tryCatch({
  # Retrieve full sample metadata
  metadata <- read.csv("~/HGM/full_metadata_prefixed_joined_simple_join.csv")
  
  # Retrieve all LCMS data
  all_merged_df <- read.csv("~/HGM/ALL_MERGED_SAMPLE_ANNOTATION.csv") %>%
    filter(confidence.level <= 3)
  
  # Retrieve all hits from immune screen
  hits <- read.csv("~/PMC266/PMC_266_CompAust_Hits.csv")
  
  CA.lib <- read.csv("~/CA/CompAus_Overlapping_CIDs.csv")
  
  Total.Annotations <- read.csv("Y:/MA_BPA_Microbiome/Total-List-Of-Annotations.csv")
  
}, error = function(e) {
  stop("Error loading one or more CSV files. Please check the file paths. Original error: ", e$message)
})

# ------------------------------------------------------------------------------#
#1 Prepare hits and all_merged_df for matching
# This block is critical for setting up the CID columns correctly before the loop.


# Perform a left join to add 'cid' from CA.lib to 'hits'.
# Ensure 'cid' column from CA.lib is named 'cid' in hits after this join.
hits <- hits %>%
  left_join(CA.lib %>% select(Compound_Name, CID), by = "Compound_Name")

# 2. Prepare all_merged_df: Add CID from Total.Annotations
# This ensures 'all_merged_df' has a 'CID' column to be used for fallback matching.
all_merged_df <- all_merged_df %>%
  left_join(Total.Annotations %>% select(compound.name, CID),
            by = "compound.name")

# Rename the 'CID' column in all_merged_df for consistency if needed.
# Your original script renames 'cid' to 'CID' but 'all_merged_df' gets 'CID' directly from Total.Annotations.
# This line might be redundant or could cause issues if 'cid' was not pre-existing in all_merged_df.
# If Total.Annotations.CID is already named 'CID', this line is not needed.
# If it was named something else (e.g., Total.Annotations.CID_Value), you'd rename it here.
# Assuming Total.Annotations$CID correctly provides 'CID', no rename needed here for 'cid' to 'CID'.
# all_merged_df <- all_merged_df %>% rename(CID = cid) # <--- Commenting out unless 'cid' truly exists and needs renaming to 'CID'

# Now, ensure that 'hits' has a 'cid' column (from CA.lib join) and 'all_merged_df' has a 'CID' column (from Total.Annotations join)

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
  left_join(metadata, by = c("Samples" = "A.HGMA.ID"))


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
for (current_hit_id in unique(hits$HitID)) {
  
  # Extract current hit row
  current_hit <- hits %>% filter(HitID == current_hit_id)
  target_compound_name <- current_hit$Compound_Name
  target_cid_from_hit <- current_hit %>% pull(CID)
  
  cat("\nProcessing Hit ID:", current_hit_id, "- Compound:", target_compound_name, "\n")
  
  # --- Try Primary Match by Compound_Name ---
  target_rows_by_compound <- all_merged_df_with_metadata %>%
    filter(!is.na(area), area > 0) %>%
    filter(compound.name == target_compound_name)
  
  if (nrow(target_rows_by_compound) > 0) {
    target_rows <- target_rows_by_compound
    match_type <- "by_Compound_Name"
    cat("  -> Found", nrow(target_rows), "rows by Compound_Name.\n")
    
  } else if (length(target_cid_from_hit) == 1 && !is.na(target_cid_from_hit)) {
    # --- Fallback Match by CID ---
    target_rows_by_cid <- all_merged_df_with_metadata %>%
      filter(!is.na(area), area > 0) %>%
      filter(CID == target_cid_from_hit)
    
    if (nrow(target_rows_by_cid) > 0) {
      target_rows <- target_rows_by_cid
      match_type <- "by_CID_fallback"
      cat("  -> No match by Compound_Name. Found", nrow(target_rows), "rows by CID (fallback).\n")
    } else {
      target_rows <- data.frame()
      match_type <- "none"
      cat("  -> No match found by Compound_Name or CID.\n")
    }
  } else {
    target_rows <- data.frame()
    match_type <- "none"
    cat("  -> No CID available for fallback and no compound match.\n")
  }
  
  # --- Export if matches were found ---
  if (nrow(target_rows) > 0) {
    # Create a safe filename
    safe_filename <- gsub("[^a-zA-Z0-9_]", "_", paste0(target_compound_name, "_HitID_", current_hit_id, "_", match_type))
    output_path <- file.path(output_dir, paste0("compound_query_", safe_filename, ".csv"))
    
    # Write CSV
    write_csv(target_rows, output_path)
    cat("  -> Data saved to:", output_path, "\n")
    
  } else {
    cat("  -> No data to export for Hit ID", current_hit_id, "\n")
  }
}

cat("Processing complete.\n")