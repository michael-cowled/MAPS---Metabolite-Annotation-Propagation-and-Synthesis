library(dplyr)
library(readr)
library(tools)
library(purrr)  # For `keep()`

# --- Initial Setup ---
data_folder <- "Y:/MA_BPA_Microbiome/Dataset-Annotations"
file_list <- list.files(data_folder, full.names = TRUE, pattern = "\\.csv$")
dataset.info <- read.csv("HGM/D - Dataset.csv", stringsAsFactors = FALSE)

# Filter to only include datasets with Quality == "Good"
good_datasets <- dataset.info %>%
  filter(Quality == "Good") %>%
  pull(HGMD.ID) %>%
  unique()

# Filter file_list to only include 'Good' datasets
file_list <- file_list %>%
  keep(~ tools::file_path_sans_ext(basename(.x)) %in% good_datasets)

all_unique_entries <- tibble(
  feature.ID = character(),
  rt = numeric(),
  mz = numeric(),
  Best.Annotation = character(),
  Best.Annotation.Smiles = character(),
  Best.Annotation.Confidence.Score = numeric(),
  Best.Annotation.Confidence.Level = numeric(),
  Best.Annotation.Type = character(),
  dataset.ID = character()
)

# --- Process Each Dataset ---
for (file_path in file_list) {
  file_name <- basename(file_path)
  dataset_id <- tools::file_path_sans_ext(file_name)
  
  dataset <- tryCatch({
    readr::read_csv(file_path, show_col_types = FALSE)
  }, error = function(e) {
    message(paste("Error reading file:", file_path, "-", e$message))
    return(NULL)
  })
  
  if (is.null(dataset)) next
  
  if (!"Best.Annotation.Smiles" %in% names(dataset)) {
    message(paste("Warning: 'Best.Annotation.Smiles' not found in", file_name, "- skipping"))
    next
  }
  
  if (!"feature.ID" %in% names(dataset)) {
    message(paste("Warning: 'feature.ID' not found in", file_name, "- skipping"))
    next
  }
  
  dataset$feature.ID <- as.character(dataset$feature.ID)
  dataset$dataset.ID <- dataset_id
  
  current_unique_entries <- dataset %>%
    filter(!is.na(Best.Annotation.Smiles), Best.Annotation.Smiles != "N/A") %>%
    filter(Best.Annotation.Confidence.Level <= 2) %>%
    distinct(Best.Annotation.Smiles, .keep_all = TRUE) %>%
    select(feature.ID, rt, mz, Best.Annotation, Best.Annotation.Smiles,
           Best.Annotation.Confidence.Score, Best.Annotation.Confidence.Level,
           Best.Annotation.Type, dataset.ID)
  
  all_unique_entries <- bind_rows(all_unique_entries, current_unique_entries)
}

# --- Join Dataset Info ---
dataset.info$HGMD.ID <- as.character(dataset.info$HGMD.ID)
all_unique_entries$dataset.ID <- as.character(all_unique_entries$dataset.ID)

all_unique_entries <- all_unique_entries %>%
  left_join(dataset.info %>% select(HGMD.ID, column.type), by = c("dataset.ID" = "HGMD.ID"))

# Optional: Check for unmatched dataset.IDs
unmatched <- setdiff(unique(all_unique_entries$dataset.ID), dataset.info$HGMD.ID)
if (length(unmatched) > 0) {
  warning("These dataset.ID values have no matching column.type: ", paste(unmatched, collapse = ", "))
}

# --- Calculate Annotation Frequency ---
# Placeholder: fix_compound_names() must be defined somewhere
final_combined_unique_entries_filtered <- fix_compound_names(all_unique_entries, "Best.Annotation")

final_data_with_frequency <- final_combined_unique_entries_filtered %>%
  group_by(Best.Annotation) %>%
  mutate(Frequency = n_distinct(dataset.ID)) %>%
  ungroup()

final_combined_unique_entries_filtered <- final_data_with_frequency %>%
  group_by(column.type) %>%
  distinct(Best.Annotation.Smiles, .keep_all = TRUE) %>%
  distinct(Best.Annotation, .keep_all = TRUE) %>%
  ungroup()

# --- Write Final CSV ---
write.csv(final_combined_unique_entries_filtered,
          "Y:/MA_BPA_Microbiome/Total-List-Of-Annotations.csv",
          row.names = FALSE)

# --- Post-standardisation with PubChem ---
final_data_after_initial_annotation_and_frequency <- final_combined_unique_entries_filtered %>%
  group_by(Best.Annotation) %>%
  mutate(Initial_Frequency = n_distinct(dataset.ID)) %>%
  ungroup()

final_combined_unique_entries_filtered <- final_data_after_initial_annotation_and_frequency %>%
  group_by(column.type) %>%
  distinct(Best.Annotation.Smiles, .keep_all = TRUE) %>%
  distinct(Best.Annotation, .keep_all = TRUE) %>%
  ungroup()

final_output_with_final_frequencies <- final_combined_unique_entries_filtered %>%
  group_by(Best.Annotation) %>%
  mutate(Final_Frequency = n_distinct(dataset.ID)) %>%
  ungroup()

write.csv(final_output_with_final_frequencies,
          "Y:/MA_BPA_Microbiome/Total-List-Of-Annotations-post-standardisation.csv",
          row.names = FALSE)
