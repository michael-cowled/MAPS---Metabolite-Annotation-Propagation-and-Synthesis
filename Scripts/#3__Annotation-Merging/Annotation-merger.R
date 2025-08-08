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

# Prepare an empty tibble to collect results
all_unique_entries <- tibble(
  feature.ID = character(),
  rt = numeric(),
  mz = numeric(),
  compound.name = character(),
  smiles = character(),
  CID = numeric(),
  confidence.score = numeric(),
  confidence.level = numeric(),
  id.prob = character(),
  annotation.type = character(),
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
  
  dataset <- dataset %>%
    mutate(
      rt = as.numeric(rt),
      mz = as.numeric(mz),
      CID = as.numeric(CID),
      confidence.score = as.numeric(confidence.score),
      confidence.level = as.numeric(confidence.level),
      feature.ID = as.character(feature.ID),
      id.prob = as.character(id.prob),
      smiles = as.character(smiles),
      compound.name = as.character(compound.name),
      annotation.type = as.character(annotation.type)
    )
  
  
  if (is.null(dataset)) next
  
  if (!"smiles" %in% names(dataset)) {
    message(paste("Warning: 'smiles' not found in", file_name, "- skipping"))
    next
  }
  
  if (!"feature.ID" %in% names(dataset)) {
    message(paste("Warning: 'feature.ID' not found in", file_name, "- skipping"))
    next
  }
  
  dataset$feature.ID <- as.character(dataset$feature.ID)
  dataset$dataset.ID <- dataset_id
  dataset$id.prob <- as.character(dataset$id.prob)  # <-- Fix type mismatch
  
  current_unique_entries <- dataset %>%
    filter(confidence.level <= 3) %>%
    filter(!grepl("analogue|candidate", compound.name, ignore.case = TRUE)) %>%
    select(feature.ID, rt, mz, compound.name, smiles, CID,
           confidence.score, id.prob, confidence.level,
           annotation.type, dataset.ID)
  
  all_unique_entries <- bind_rows(all_unique_entries, current_unique_entries)
}


# --- Join Dataset Info ---
dataset.info$HGMD.ID <- as.character(dataset.info$HGMD.ID)
all_unique_entries$dataset.ID <- as.character(all_unique_entries$dataset.ID)

all_unique_entries <- all_unique_entries %>%
  left_join(dataset.info %>% select(HGMD.ID, column.type), by = c("dataset.ID" = "HGMD.ID"))

# Remove PUBCHEM fallback labels
all_unique_entries <- all_unique_entries %>%
  filter(!grepl("PUBCHEM", compound.name, ignore.case = TRUE))

# Optional: Check for unmatched dataset.IDs
unmatched <- setdiff(unique(all_unique_entries$dataset.ID), dataset.info$HGMD.ID)
if (length(unmatched) > 0) {
  warning("These dataset.ID values have no matching column.type: ", paste(unmatched, collapse = ", "))
}

# --- Calculate Annotation Frequency ---
# Placeholder: fix_compound_names() must be defined elsewhere

final_combined_unique_entries_filtered <- fix_compound_names(all_unique_entries, "compound.name")

final_data_with_frequency <- final_combined_unique_entries_filtered %>%
  group_by(compound.name) %>%
  mutate(Frequency = n_distinct(dataset.ID)) %>%
  ungroup()

compound_frequencies <- final_data_with_frequency %>%
  filter(!is.na(compound.name)) %>%
  group_by(compound.name) %>%
  summarise(Frequency = n_distinct(dataset.ID))

# Step 1: Assign column-type priority
priority_filtered <- final_data_with_frequency %>%
  mutate(priority = case_when(
    column.type %in% c("Phe-Hex", "HILIC") ~ 1,
    column.type == "C18" ~ 2,
    TRUE ~ 3
  ))

# Step 2: Keep only best annotation per (smiles + column.type)
best_per_smiles_coltype <- priority_filtered %>%
  arrange(confidence.level, desc(id.prob), desc(confidence.score)) %>%
  group_by(compound.name, column.type) %>%
  slice_head(n = 1) %>%
  ungroup()

# Step 3: Drop C18 if Phe-Hex exists for that compound
final_combined_unique_entries_filtered <- best_per_smiles_coltype %>%
  group_by(compound.name) %>%
  filter(
    column.type %in% c("Phe-Hex", "HILIC") |
      (column.type == "C18" & !any(column.type == "Phe-Hex"))
  ) %>%
  ungroup()

# Step 4: Add frequency of datasets in which each compound was detected
final_combined_unique_entries_filtered <- final_combined_unique_entries_filtered %>%
  group_by(compound.name) %>%
  mutate(Frequency = n_distinct(dataset.ID)) %>%
  ungroup()

# Step 5: Join the true frequency back to the filtered results
final_combined_unique_entries_filtered <- final_combined_unique_entries_filtered %>%
  left_join(compound_frequencies, by = "compound.name")

#Final tidying

final_combined_unique_entries_filtered <- final_combined_unique_entries_filtered %>%
  select(-Frequency.x, -priority) %>%
  filter(!is.na(compound.name), !is.na(smiles))

# --- Write Final CSV ---
write_csv(final_combined_unique_entries_filtered,
          "Y:/MA_BPA_Microbiome/Total-List-Of-Annotations-v3.csv",
          row.names = FALSE)

hilic <- filter(final_combined_unique_entries_filtered, column.type == "HILIC")
phehex <- filter(final_combined_unique_entries_filtered, column.type == "Phe-Hex" | column.type == "C18")

write_csv(hilic,
          "Y:/MA_BPA_Microbiome/HILIC-List-Of-Annotations.csv",
          row.names = FALSE)

write_csv(phehex,
          "Y:/MA_BPA_Microbiome/PheHex-List-Of-Annotations.csv",
          row.names = FALSE)

# --- USI-style Counts ---
usi_count_func <- function(dataset) {
  dataset %>%
    distinct(compound.name, confidence.level) %>%  # No need for .keep_all = TRUE
    group_by(confidence.level) %>%
    summarize(total_usi = n(), .groups = "drop")
}

cat("USI counts (Total):\n")
print(usi_count_func(final_combined_unique_entries_filtered))

cat("USI counts (HILIC):\n")
print(usi_count_func(hilic))

cat("USI counts (Phe-Hex/C18):\n")
print(usi_count_func(phehex))
