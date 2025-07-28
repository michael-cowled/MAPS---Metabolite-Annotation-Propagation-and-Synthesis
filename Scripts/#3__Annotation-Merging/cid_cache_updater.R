# Load the dplyr library for data manipulation
# If you don't have it, run: install.packages("dplyr")
library(dplyr)

# --- 1. SETUP ---
# 📂 Replace these with the actual file paths for your two files
path_to_file1 <- "~/cid_cache.csv"
path_to_file2 <- "C:/Users/mcowled/Downloads/cid_cache (2).csv"

# --- 2. READ DATA ---
# Read both CSV files into R data frames
# The tryCatch block will stop the script if a file doesn't exist
file1_data <- tryCatch({
  read.csv(path_to_file1)
}, error = function(e) {
  stop("Error: Could not read the first file. Please check the path: ", path_to_file1)
})

file2_data <- tryCatch({
  read.csv(path_to_file2)
}, error = function(e) {
  stop("Error: Could not read the second file. Please check the path: ", path_to_file2)
})

# --- 3. FIND NEW ROWS ---
# Use an anti_join to find rows that are in the second file but NOT in the first.
# This compares entire rows to find what's new.
new_rows <- anti_join(file2_data, file1_data)

# --- 4. APPEND AND SAVE ---
# Check if there are actually any new rows to add
if (nrow(new_rows) > 0) {
  
  # ✨ Announce how many new rows were found
  cat("Found", nrow(new_rows), "new rows to add.\n")
  
  # Combine the original data from the first file with the new rows
  updated_file1_data <- bind_rows(file1_data, new_rows)
  
  # 💾 Save the updated data frame back to the first file's location
  # This will overwrite the original file with the newly combined data.
  # row.names = FALSE is crucial to prevent R from adding an extra column.
  write.csv(updated_file1_data, path_to_file1, row.names = FALSE)
  
  cat("Success! The first file has been updated.\n")
  
} else {
  
  # Announce that no differences were found
  cat("No new rows found in the second file. The first file is already up-to-date.\n")
  
}