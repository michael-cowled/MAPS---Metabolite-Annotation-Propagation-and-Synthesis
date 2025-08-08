######################### USER GUIDE #########################################

#### This script is to extract all the sample df and annotation df in the LCMS folder in Y drive to bind all information together   ####

#### Tip: could exclude some merged_df from certain folders in line 95 ###############


library(dplyr)
library(readr)
library(purrr)
library(tidyr)
library(readr)
library(stringr)

# Retrieve full metadata
metadata<-read.csv("~/HGM/full_metadata_prefixed_joined_simple_join.csv")

# ---------------------Retrieve sample files and annotation files from mediaflux -----------------------------------####

# Set paths to the two central folders
annotation_path <- "Y:/MA_BPA_Microbiome/Dataset-Annotations"
sample_path <- "Y:/MA_BPA_Microbiome/Dataset-Abundances"

# List annotation and sample files
annotation_files <- list.files(annotation_path, pattern = "^HGMD_\\d{4}\\.csv$", full.names = TRUE)
sample_files <- list.files(sample_path, pattern = "^HGMD_\\d{4}-samples-df\\.csv$", full.names = TRUE)

# Extract HGMD IDs from filenames (without path)
annotation_ids <- str_extract(basename(annotation_files), "HGMD_\\d{4}")
sample_ids <- str_extract(basename(sample_files), "HGMD_\\d{4}")

# Get the intersection of IDs (i.e., those that have both annotation and sample files)
common_ids <- intersect(annotation_ids, sample_ids)

# Loop through each valid HGMD_ID
for (hgmd_id in common_ids) {
  
  # Build file paths
  annotation_file <- file.path(annotation_path, paste0(hgmd_id, ".csv"))
  sample_file <- file.path(sample_path, paste0(hgmd_id, "-samples-df.csv"))
  
  # Read and assign data frames
  assign(paste0("df_", hgmd_id, "_annotation"), read_csv(annotation_file, show_col_types = FALSE))
  assign(paste0("df_", hgmd_id, "_samples"), read_csv(sample_file, show_col_types = FALSE))
}

# Optional: list all loaded objects
ls(pattern = "^df_HGMD_\\d{4}_")
###-----------------------------------Another loop to merge all pairs-----------------------------------------------------##########

# Get all sample dataframe names
sample_dfs <- ls(pattern = "^df_HGMD_\\d{4}_samples$")

for (sample_name in sample_dfs) {
  hgmd_id <- str_extract(sample_name, "HGMD_\\d{4}")
  annotation_name <- paste0("df_", hgmd_id, "_annotation")
  
  if (exists(sample_name) && exists(annotation_name)) {
    tryCatch({
      # Get sample df and rename 'samples' to 'Samples'
      df_sample <- get(sample_name) %>% 
        rename(Samples = samples)
      
      # Get annotation df and remove conflicting columns and index column
      df_annotation <- get(annotation_name) %>% 
        select(-feature.usi, -any_of(c("Samples", "...1")))
      
      message("Merging HGMD: ", hgmd_id)
      
      # Join dataframes on feature.ID
      merged_df <- df_sample %>% 
        left_join(df_annotation, by = "feature.ID") %>%
        
        # Resolve duplicated compound.name and smiles columns
        mutate(
          compound.name = coalesce(compound.name.y, compound.name.x),
          smiles = coalesce(smiles.y, smiles.x)
        ) %>%
        # Drop original duplicated columns
        select(-compound.name.x, -compound.name.y, -smiles.x, -smiles.y) %>%
        
        # Select relevant columns with safety using any_of()
        select(any_of(c(
          "feature.ID",
          "Samples",
          "feature.usi",
          "area",
          "rt",
          "mz",
          "compound.name",
          "smiles",
          "confidence.score",
          "NPC.superclass",
          "confidence.level",
          "id.prob",
          "annotation.type"
        ))) %>%
        
        mutate(HGMD.ID = hgmd_id)
      
      # Remove redundant index column if present
      if (colnames(merged_df)[1] == "...1") {
        merged_df <- merged_df[, -1, drop = FALSE]
      }
      
      # Assign merged dataframe to global env with a new name
      assign(paste0("merged_df_", hgmd_id), merged_df, envir = .GlobalEnv)
      
    }, error = function(e) {
      message(paste("❌ Unable to merge", hgmd_id, ":", e$message))
    })
  }
}

###Bind all merged df#####
##############---------------------------------------------------------------------------------------------------##############
# Get all merged data frame names
merged_dfs <- ls(pattern = "^merged_df_HGMD_\\d{4}$")

# Filter out those in the exclude list
merged_dfs <- merged_dfs[!str_extract(merged_dfs, "HGMD_\\d{4}") %in% exclude_ids]

# Bind all into one data frame
all_merged_df <- bind_rows(
  lapply(merged_dfs, get),
  .id = "source_df"
)

# Check if (Samples, HGMD.ID, feature.usi) pairs are unique
is_unique <- all_merged_df %>%
  count(Samples, HGMD.ID, feature.usi) %>%
  filter(n > 1)

# View duplicates (if any)
if (nrow(is_unique) > 0) {
  print(is_unique)
} else {
  message("No duplicate (Samples, HGMD.ID, feature.usi) combinations found.")
}

write_csv(all_merged_df, "~/HGM/ALL_MERGED_SAMPLE_ANNOTATION.csv")