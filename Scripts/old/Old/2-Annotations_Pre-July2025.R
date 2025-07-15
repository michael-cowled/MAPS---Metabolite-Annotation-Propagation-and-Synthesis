#### ANNOTATION TABLE SCRIPT FOR AUSTRALIAN HUMAN GUT METABOLOME DATABASE ####

#!# Run 1-Functions first! #!!#

#------------------------------------------------------------------------------#
              ###---START OF USER-FED INFORMATION---###

# Dataset ID: From Data Management Plan:

  dataset.id <- "HGMD_0108"     #####Change to dataset of interest#####

# Specify the path to the Data Management Plan (unless moved to Mediaflux then should be consistent)

  excel_file <- "C:/Users/mcowled/The University of Melbourne/MA Human Gut Metabolome - Documents/Data Management and Analysis/HGM - Data Management System.xlsx"

# Annotation Acceptance Probabilities (Defaults for Exploratory Analyses)
# Increase stringency if the application demands it.

  gnps.prob <- 0.7              # Default is 0.7 (and should be changed to those set in run)
  canopus.prob <- 0.7           # Default is 0.7
  csi.prob <- 0.64            # Default is 0.64 ##10% FDR
  ms2query.prob <- 0.63         # Default is 0.63 ##Recommended by paper
  
# Specify rt tolerance of standards (min)
  
  rt.tol <- 0.1               # Default is 0.1 min for C18 (use 0.2 min for HILIC)

              ###---END OF USER-FED INFORMATION---###
#------------------------------------------------------------------------------#
  
  
## 1. Loading libraries and Cache
required_packages <- c("dplyr", "tidyr", "stringr", "readr", "reshape2", "ggplot2", 
                         "svglite", "readxl", "openxlsx", "tidyverse", "rvest", "jsonlite", "xml2")
check_and_install(required_packages)

# --- Load or Initialize Cache ---
cid_cache_file <- "cid_cache.csv"
if (file.exists(cid_cache_file)) {
  cid_cache_df <- read.csv(cid_cache_file, stringsAsFactors = FALSE)
} else {
  cid_cache_df <- data.frame(
    LookupName = character(),
    ResolvedName = character(),
    SMILES = character(),
    CID = numeric(),
    MolecularFormula = character(),
    MonoisotopicMass = numeric(),
    stringsAsFactors = FALSE
  )
}

## 2. Updating Metadata and extract gnps task ID                                ----> Remove from published version
sheet_names <- excel_sheets(excel_file) # Get the sheet names
for (sheet in sheet_names) {
  data <- read_excel(excel_file, sheet = sheet)
  write.csv(data, file = paste0("HGM/", sheet, ".csv"), row.names = FALSE)
}

dataset <- read.csv("HGM/D - Dataset.csv") %>%
  filter(HGMD.ID == dataset.id)
gnps.task.id <- dataset$gnps.task.ID[1] #                                       ----> Put gnps.task.ID as a required user-fed info in published

if (is.na(gnps.task.id)) {
  stop("gnps.task.ID is missing. Add to 'dataset' csv before re-running.")  # Stop execution
}

## 3. File List Extractor                                                       ----> Remove from published version
dataset <- read.csv("HGM/D - Dataset.csv") %>%
  filter(HGMD.ID == dataset.id)
output_directory <- dataset$Processed.Data.Folder[1] # Take the first value if there are multiple
folder <- paste0(dataset$Processed.Data.Folder, "\\mzml")
print(folder)
if (dir.exists(folder)) {
  files <- list.files(folder)
  full_paths <- file.path(folder, files)
  full_paths_windows <- gsub("/", "\\\\", full_paths)
  file_info <- data.frame(
    folder = folder,
    filename = files,
    full_path = full_paths_windows
  )
  print(file_info)
  output_filename <- paste0(dataset.id, "_file_list.csv")
  output_path <- file.path(output_directory, output_filename) # changed here
  write.csv(file_info, file = output_path, row.names = FALSE) # and here
  print(paste("File list saved to:", output_path)) # and here
} else {
  print("Folder does not exist. Check the path.")
}

## 4. Processed Data Check
#Check if all data is present in folder with correct naming conventions         ###CHECK SPELLING!!!!###
#Folder with mzmine, ms2query, sirius and gnps results
folder <- dataset$Processed.Data.Folder[1]  ##Determines folder to process based on dataset.id ---->Put as required user-input information

# Construct the full paths for each data file
mzmine.data <- paste0(folder, "/mzmine/ms1-and-ms2.csv")
mzmine.annotations <- paste0(folder, "/mzmine/data_annotations.csv")
canopus.data <- paste0(folder, "/sirius/canopus_structure_summary.tsv")
csi.data <- paste0(folder, "/sirius/structure_identifications_top-100.tsv")
zodiac.data <- paste0(folder, "/sirius/formula_identifications.tsv")
ms2query.data <- paste0(folder, "/ms2query/ms2query.csv")
cytoscape <- paste0(folder, "/gnps/cytoscape.csv")

validate_and_get_paths(folder)

## 5. Load in MZMINE data
mzmine.annotations <- read.csv(mzmine.annotations) %>%
  # First, ensure distinct compound_name per id
  distinct(id, compound_name, .keep_all = TRUE) %>%
  group_by(id) %>%
  # Remove all rows with the lowest score per id
  filter(score > min(score) | n() == 1) %>%
  ungroup() %>%
  # Now, reduce to only one row per compound_name (highest score overall)
  group_by(compound_name) %>%
  filter(score == max(score)) %>%
  slice(1) %>%  # If there's a tie on score, pick the first row arbitrarily
  ungroup()

#Calculating ID probability of level 1 annotations
mzmine.annotations$mzmine.ID.prob <- NA
#Standardisation of compound names (retrieving from a locally stored cache/pubchem)
# --- Main Processing Loop ---
mzmine.annotations <- standardise_annotation(mzmine.annotations, "compound_name", "smiles", "authentic.standard")
write.csv(cid_cache_df, cid_cache_file, row.names = FALSE)

# --- Populate from cache ---
for (i in seq_len(nrow(mzmine.annotations))) {
  current_name <- mzmine.annotations$compound_name[i]
  row_match <- cid_cache_df[cid_cache_df$LookupName == current_name, ]
  
  if (nrow(row_match) > 0) {
    mzmine.annotations$authentic.standard.CID[i] <- row_match$CID[1]
    mzmine.annotations$smiles[i] <- row_match$SMILES[1]
    mzmine.annotations$mol_formula[i] <- row_match$MolecularFormula[1]
    mzmine.annotations$authentic.standard.MonoisotopicMass[i] <- row_match$MonoisotopicMass[1]
  }
}

mzmine.annotations.final <- mzmine.annotations %>%
  group_by(id) %>%
  # Add mzmine.ID.prob = 1 / number of rows per id
  mutate(mzmine.ID.prob = 1 / n()) %>%
  # Keep only the row(s) with the highest score per id
  filter(score == max(score)) %>%
  slice(1) %>%  # In case of ties, keep one arbitrarily
  ungroup() %>%
  # Keep only relevant columns:
  select(id, compound_name, score, smiles, mzmine.ID.prob, mol_formula, authentic.standard.CID, authentic.standard.MonoisotopicMass)
names(mzmine.annotations.final) <- c('feature.ID', "authentic.standard", "authentic.standard.score", 
                                     "authentic.standard.smiles", "authentic.standard.ID.prob",
                                     "authentic.standard.molecular.formula", "authentic.standard.CID", "authentic.standard.MonoisotopicMass")
mzmine.annotations.final$feature.ID <- as.numeric(mzmine.annotations.final$feature.ID)
mzmine.annotations.final$confidence.level <- "1"

#Now appending to all mzmine features
mzmine.data <- read.csv(mzmine.data) # Derived from Export to CSV file (modular)
sample.data <- mzmine.data # A copy to be used for different processing later
if (!("spectral_db_matches.compound_name" %in% names(mzmine.data))) {
  mzmine.data$spectral_db_matches.compound_name <- NA
}
mzmine.data <- select(mzmine.data, "id", "rt", "mz", "ion_identities.iin_id") %>%
  mutate(across(everything(), as.character)) %>%
  mutate(across(everything(), ~ na_if(., "")))
names(mzmine.data) <- c('feature.ID', "rt", "mz", "ion.identity.ID")
mzmine.data$feature.ID <- as.numeric(mzmine.data$feature.ID)
mzmine.data <- mzmine.data %>%
  left_join(mzmine.annotations.final, by = "feature.ID")

## 6. Load, tidy and standardise GNPS2 data - works for v0.1.2, no metadata required
gnps.annotation.data <- read_tsv(paste0("https://gnps2.org/resultfile?task=", gnps.task.id, "&file=nf_output/library/merged_results_with_gnps.tsv"))
gnps.annotation.data <- gnps.annotation.data[, c(2, 4, 5, 8, 9, 15, 27, 35, 43, 45, 46)]
names(gnps.annotation.data) <- c("feature.ID", "gnps.library.name", "gnps.cosine.score", 
                                 "gnps.diff.ppm", "gnps.shared.peaks", "gnps.compound.name", 
                                 "gnps.smiles", "gnps.library.quality", "gnps.NPC.superclass", "gnps.NPC.pathway", "library.usi")
gnps.annotation.data$library.usi <- str_replace(
  gnps.annotation.data$library.usi,
  "mzspec:GNPS:GNPS-LIBRARY:(.*)",
  "mzspec:GNPS:GNPS-LIBRARY:accession:\\1"
)
gnps.annotation.data <- gnps.annotation.data[, c(1, 6, 7, 3:5, 2, 8, 10, 9, 11)] #Reorders columns
gnps.cluster.data <- read_tsv(paste0("https://gnps2.org/resultfile?task=", gnps.task.id, "&file=nf_output/networking/clustersummary_with_network.tsv"))
gnps.cluster.data <- select(gnps.cluster.data,  'cluster index', 'component')
names(gnps.cluster.data) <- c('feature.ID', "gnps.cluster.ID")
gnps.data <- gnps.cluster.data %>%
  full_join(gnps.annotation.data, by = "feature.ID") %>%
  filter(gnps.diff.ppm <= 5) %>%
  filter(gnps.cosine.score >= gnps.prob)
gnps.data$feature.usi <- paste0("mzspec:GNPS2:TASK-", gnps.task.id, "-nf_output/clustering/spectra_reformatted.mgf:scan:", gnps.data$feature.ID) #Adds unique spectra identifiers
gnps.cluster.pairs <- read_tsv(paste0("https://gnps2.org/resultfile?task=", gnps.task.id, "&file=nf_output/networking/filtered_pairs.tsv"))
gnps.data$gnps.in.silico.bile.acid.info <- NA
gnps.data$gnps.in.silico.bile.acid.info[grepl("Candidate ", gnps.data$gnps.compound.name)] <- 
  gnps.data$gnps.compound.name[grepl("Candidate ", gnps.data$gnps.compound.name)]
cols_needed <- c("gnps.CID", "gnps.smiles", "gnps.MolecularFormula", "gnps.MonoisotopicMass")
for (col in cols_needed) {
  if (!col %in% names(gnps.data)) {
    gnps.data[[col]] <- if (col == "gnps.MonoisotopicMass") NA_real_ else NA_character_
  }
}
# Ensure required output columns exist
gnps.data[[paste0("gnps.CID")]] <- NA
gnps.data[[paste0("gnps.MolecularFormula")]] <- NA
gnps.data[[paste0("gnps.MonoisotopicMass")]] <- NA

#Create level 2, and tidy
gnps.data.lv2 <- filter(gnps.data.lv2, gnps.library.quality != "Insilico")
gnps.data.lv2$confidence.level <- "2"
gnps.data.lv2 <- fix_compound_names(gnps.data.lv2, "gnps.compound.name")

#Create level 3, and tidy
gnps.data.lv3 <- filter(gnps.data, gnps.library.quality == "Insilico")
gnps.data.lv3$confidence.level <- "3"
gnps.data.lv3 <- gnps.data.lv3 %>%
  mutate(
    gnps.compound.name = ifelse(
      grepl("Candidate ", gnps.compound.name),
      sub("\\s*\\(.*$", "", gnps.compound.name),
      gnps.compound.name
    )
  )
gnps.data.lv3$gnps.in.silico.bile.acid.info <- gsub('"', "", gnps.data.lv3$gnps.in.silico.bile.acid.info)     
gnps.data.lv3$gnps.compound.name <- gsub('"', "", gnps.data.lv3$gnps.compound.name)       

#Standardisation of compound names (retrieving from a locally stored cache/pubchem)
gnps.data.lv2 <- deduplicate_data(gnps.data.lv2, gnps.compound.name, gnps.cosine.score) #Pre-standardisation filtering for duplicates
gnps.data.lv2 <- standardise_annotation(gnps.data.lv2, "gnps.compound.name", "gnps.smiles", "gnps") #Standardisation of naming using pubchem, remove for time-savings
gnps.data.lv2 <- deduplicate_data(gnps.data.lv2, , gnps.compound.name, gnps.cosine.score) #Post-standardisation filtering of duplicates
write.csv(cid_cache_df, cid_cache_file, row.names = FALSE)

## 7. Load in MS2QUERY data
ms2query.data <- read.csv(ms2query.data)
ms2query.data  <- ms2query.data[, c(2, 3, 7, 8, 10, 17, 18)]
names(ms2query.data) <- c("ms2query.score", "ms2query.mzdiff", 
                          "ms2query.analogue.compound.name", "ms2query.smiles", 
                          'feature.ID', "ms2query.NPC.superclass", "ms2query.NPC.pathway")
ms2query.data <- ms2query.data[, c(3, 1, 2, 4, 5, 7, 6)] #Reorders columns
ms2query.data <- fix_compound_names(ms2query.data, "ms2query.analogue.compound.name")  ##df, and column to fix

#Create level 2, and tidy
ms2query.data.lv2 <-ms2query.data %>%
  filter(ms2query.score > ms2query.prob) %>%
  filter(ms2query.mzdiff <= 0.001) 
ms2query.data.lv2$confidence.level <- "2"
ms2query.data.lv2 <- deduplicate_data(ms2query.data.lv2, ms2query.analogue.compound.name, ms2query.score) #Pre-standardisation filtering for duplicates
ms2query.data.lv2 <- standardise_annotation(ms2query.data.lv2, "ms2query.analogue.compound.name", "ms2query.smiles", "ms2query") #Standardisation of naming using pubchem, remove for time-savings
ms2query.data.lv2  <- deduplicate_data(ms2query.data.lv2, ms2query.analogue.compound.name, ms2query.score) #Post-standardisation filtering of duplicates
write.csv(cid_cache_df, cid_cache_file, row.names = FALSE)

#Create level 3, and tidy
ms2query.data.lv3 <- ms2query.data %>%
  filter(ms2query.score > ms2query.prob) %>%
  filter(ms2query.mzdiff >= 0.001) %>%
  filter(!feature.ID %in% mzmine.annotations.final$feature.ID) #removes any annotations that conflict with level 1 annotations
ms2query.data.lv3$confidence.level <- "3"

## 8. Level 2 ID probabilities
gnps.data <- filter_data_by_standards(gnps.data, mzmine.data, "gnps.compound.name", "feature.ID", rt.tol)
ms2query.data.lv2 <- filter_data_by_standards(ms2query.data.lv2, mzmine.data, "ms2query.analogue.compound.name", "feature.ID", rt.tol)

#-----------------#
###AFTER LUNCH
gnps.data <- gnps.data %>%
  mutate(gnps.ID.prob = 0)
gnps.data <- gnps.data %>%
  group_by(feature.ID) %>%
  mutate(
    n_above_thresh = sum(gnps.cosine.score >= gnps.prob),
    gnps.ID.prob = ifelse(gnps.cosine.score >= gnps.prob & n_above_thresh > 0, 1 / n_above_thresh, 0)
  ) %>%
  ungroup() %>%
  select(-n_above_thresh)
gnps.data <- gnps.data %>%
  group_by(feature.ID) %>%
  arrange(desc(gnps.cosine.score)) %>%
  slice(1) %>%
  ungroup()    
    
## 9. Load in SIRIUS data:      Compatible with v6.1.0 onwards
canopus.data <- read_tsv(canopus.data)
canopus.data <- canopus.data[, c(5:8, 27)]
names(canopus.data) <- c("canopus.NPC.pathway", "canopus.NPC.pathway.probability", 
                         "canopus.NPC.superclass", "canopus.NPC.superclass.probability", 
                         'feature.ID')
canopus.data <- canopus.data %>%
  group_by(feature.ID) %>%
  filter(!(all(canopus.NPC.pathway.probability == 0))) %>%  # Remove groups where all scores are 0
  filter(canopus.NPC.pathway.probability == max(canopus.NPC.pathway.probability, na.rm = TRUE)) %>%  # Keep the highest score
  slice(1) %>%  # Arbitrarily keep one if multiple rows have the same non-zero score
  ungroup()

csi.data <- read_tsv(csi.data)
csi.data <- csi.data[, c(3, 14, 15, 25)]
names(csi.data) <- c("csi.confidence.score", "csi.compound.name", "csi.smiles", 'feature.ID')
csi.data <- csi.data[, c(2, 1, 3:ncol(csi.data))] #Swap cols 1 and 2
# Initial fixes
csi.data$csi.confidence.score[csi.data$csi.confidence.score == -Inf] <- 0
csi.data$csi.compound.name[grepl("Solaparnaine", csi.data$csi.compound.name, ignore.case = TRUE)] <- "Solaparnaine" ##troublesome case
csi.data$csi.confidence.score <- as.numeric(csi.data$csi.confidence.score)

# Calculation of Identification Probability (Wishart paper)
csi.data <- csi.data %>%
  mutate(csi.ID.prob = 0)
csi.data <- csi.data %>%
  group_by(feature.ID) %>%
  mutate(
    n_above_thresh = sum(csi.confidence.score >= csi.prob),
    csi.ID.prob = ifelse(csi.confidence.score >= csi.prob & n_above_thresh > 0, 1 / n_above_thresh, 0)
  ) %>%
  ungroup() %>%
  select(-n_above_thresh)
csi.data <- csi.data %>%
  group_by(feature.ID) %>%
  arrange(desc(csi.confidence.score)) %>%
  slice(1) %>%
  ungroup()  

zodiac.data <- read_tsv(zodiac.data)
zodiac.data  <- zodiac.data[, c(2, 5, 21)]
names(zodiac.data) <- c("zodiac.formula", "zodiac.score", 'feature.ID')
# Filter and keep only unique feature.ID (in the case where multiple annotations are provided)
zodiac.data <- zodiac.data %>%
  group_by(feature.ID) %>%
  filter(!(all(zodiac.score == 0))) %>%  # Remove groups where all scores are 0
  filter(zodiac.score == max(zodiac.score, na.rm = TRUE)) %>%  # Keep the highest score
  slice(1) %>%  # Arbitrarily keep one if multiple rows have the same non-zero score
  ungroup()

## 11. Merge annotations into one big table                                     
full.annotation.data <- mzmine.data %>%
  full_join(gnps.data, by = "feature.ID") %>%
  full_join(canopus.data, by = "feature.ID") %>%
  full_join(csi.data, by = "feature.ID") %>%
  full_join(zodiac.data, by = "feature.ID") %>%
  full_join(ms2query.data, by = "feature.ID")

full.annotation.data <- standardise_annotation(full.annotation.data, "authentic.standard", "authentic.standard.smiles", "authentic.standard") %>%
  select(-CID)
write.csv(cid_cache_df, cid_cache_file, row.names = FALSE)

## 12. Propagation of annotations

#Create a summary.annotation.data table, that is filtered for all annotations satisfying a certain degree of confidence
summary.annotation.data <- full.annotation.data

#Apply the conditions
summary.annotation.data$csi.compound.name[summary.annotation.data$csi.confidence.score < csi.prob] <- NA
summary.annotation.data$csi.smiles[summary.annotation.data$csi.confidence.score < csi.prob] <- NA
summary.annotation.data$ms2query.analogue.compound.name[summary.annotation.data$ms2query.score < ms2query.prob] <- NA
summary.annotation.data$ms2query.analogue.compound.name[summary.annotation.data$ms2query.mzdiff > 0.001] <- NA

# Initialize columns in propagated.annotation.data
propagated.annotation.data <- full.annotation.data %>%
  mutate(
    Probable.Analogue.Of = NA,
    Propagated.Feature.ID = NA,
    Propagated.Annotation.Type = NA, # Initialize the new column
    Propagated.Annotation.Class = NA # This is already here
  )

# Identify the rows where all specified columns contain NA
na.rows <- summary.annotation.data %>%
  filter(
    is.na(gnps.compound.name) &
      is.na(csi.compound.name) &
      is.na(csi.smiles) &
      is.na(ms2query.analogue.compound.name)
  )

# Extract the 'feature.ID' column from these rows
na.feature.ids <- na.rows$feature.ID

# Iterate over na.feature.ids to populate propagated.annotation.data
for (i in na.feature.ids) {
  paired_values <- paired_feature_finder(i)
  print(paste("Paired values for feature ID", i, ":", paste(paired_values, collapse = ", ")))
  
  # Initialize variables to store results and corresponding paired value
  result_data <- list(value = NA, column = NA, superclass = NA) # Initialize superclass here
  selected_paired_value <- NA
  
  for (value in paired_values) {
    result_data <- get_result(value, summary.annotation.data)
    if (!is.na(result_data$value)) {
      selected_paired_value <- value # Store the paired value used for the result
      break # Exit the loop if a valid result is found
    }
  }
  
  # Update the new columns in propagated.annotation.data
  propagated.annotation.data <- propagated.annotation.data %>%
    mutate(
      Probable.Analogue.Of = case_when(
        feature.ID == i ~ result_data$value,
        TRUE ~ Probable.Analogue.Of
      ),
      Propagated.Feature.ID = case_when(
        feature.ID == i ~ selected_paired_value,
        TRUE ~ Propagated.Feature.ID
      ),
      Propagated.Annotation.Type = case_when(
        feature.ID == i ~ result_data$column,
        TRUE ~ Propagated.Annotation.Type
      ),
      Propagated.Annotation.Class = case_when( # Add this new case_when
        feature.ID == i ~ result_data$superclass,
        TRUE ~ Propagated.Annotation.Class
      )
    )
}

# Calculating the mz difference of the analogue to the propagated feature:
# Convert mz column to numeric if it's not already
if (!is.numeric(summary.annotation.data$mz)) {
  summary.annotation.data$mz <- as.numeric(summary.annotation.data$mz)
}

# Add a new column for the m/z difference
propagated.annotation.data <- propagated.annotation.data %>%
  mutate(Propagated.Annotation.mz.Diff = NA)

# Iterate over the rows of propagated.annotation.data
for (i in 1:nrow(propagated.annotation.data)) {
  # Get the feature ID and the propagated feature ID
  feature_id <- propagated.annotation.data$feature.ID[i]
  propagated_feature_id <- propagated.annotation.data$Propagated.Feature.ID[i]
  
  # If a propagated feature ID exists
  if (!is.na(propagated_feature_id)) {
    # Get the m/z values for both feature IDs
    mz_value <- summary.annotation.data$mz[summary.annotation.data$feature.ID == feature_id]
    propagated_mz_value <- summary.annotation.data$mz[summary.annotation.data$feature.ID == propagated_feature_id]
    
    # Calculate the m/z difference
    mz_diff <- mz_value - propagated_mz_value
    
    # Update the new column
    propagated.annotation.data$Propagated.Annotation.mz.Diff[i] <- mz_diff
  }
}

# Inside the loop:
if (!is.na(propagated_feature_id)) {
  mz_value <- summary.annotation.data$mz[summary.annotation.data$feature.ID == feature_id]
  propagated_mz_value <- summary.annotation.data$mz[summary.annotation.data$feature.ID == propagated_feature_id]
  
  if (!is.na(mz_value) && !is.na(propagated_mz_value)) {  # Check for NA
    mz_diff <- mz_value - propagated_mz_value
    propagated.annotation.data$Propagated.Annotation.mz.Diff[i] <- mz_diff
  }
}

## 13. Sample List Data
# Now to append sample.list corresponding to each feature
# Select columns containing ".area" and rename the first column
sample.data <- sample.data %>%
  select(id, contains(".area"))

# Rename the first column to "feature.ID"
colnames(sample.data)[1] <- "feature.ID"

# Get the column names of the data frame
colnames_sample <- colnames(sample.data)

# Remove the prefix "datafile." and the suffix ".mzML.Peak.area" if present
colnames_sample <- sub("^datafile\\.", "", colnames_sample)  # Remove the prefix
colnames_sample <- sub("\\.mzML\\.area$", "", colnames_sample)  # Remove the suffix
colnames(sample.data) <- colnames_sample

sample.data2 <- sample.data ##For other uses later 

# Replace all positive areas with 1 i.e. presence.absence
sample.data[, 2:ncol(sample.data)] <- lapply(sample.data[, 2:ncol(sample.data)], function(x) {
  x[x > 0] <- 1
  return(x)
})

# Create a new column with presence list separated by semicolons
sample.data$Samples <- apply(sample.data[, 2:ncol(sample.data)], 1, function(row) {
  paste(colnames(sample.data[, 2:ncol(sample.data)])[which(row == 1)], collapse = "; ")
})

# Select only the first and the last columns (which includes the new 'Samples' column)
sample.data <- sample.data[, c(1, ncol(sample.data))]

#Append Sample data
propagated.annotation.data.with.samples <- propagated.annotation.data %>%
  full_join(sample.data, by = "feature.ID") 

## 14. Collapsing Ion Identity Networks
#performed in such a way to retain the best annotation (if multiple)

#A new editing df where features are sequentially removed:
propagation.df <- propagated.annotation.data.with.samples

#A new df with "good" annotations extracted in a similar manner to unpropagated version.
new.summary.annotation.data <- propagated.annotation.data.with.samples
#Apply the conditions
new.summary.annotation.data$csi.compound.name[new.summary.annotation.data$csi.confidence.score < csi.prob] <- NA
new.summary.annotation.data$csi.smiles[new.summary.annotation.data$csi.confidence.score < csi.prob] <- NA
new.summary.annotation.data$ms2query.analogue.compound.name[new.summary.annotation.data$ms2query.score < ms2query.prob] <- NA

iin.features <- filter(new.summary.annotation.data, !is.na(ion.identity.ID))

results <- process_all_features(data, "ion.identity.ID", c("feature1", "feature2"))

# Print results for each ID
for (id in names(results)) {
  cat("Results for ion.identity.ID =", id, ":\n")
  print(results[[id]])
}

# Update propagation.df by removing rows based on results and retaining NA IDs
propagation.df <- update_data_frame(propagation.df, results, "ion.identity.ID")

## 15. Best Annotation Extraction
#Simultaneously creates a tidied cytoscape file with the best annotation (including analogues)
#Represented as the output.

cytoscape <- read.csv(paste0(folder, "/gnps/cytoscape.csv"))

level1.annotation <- filter(propagation.df, !is.na(authentic.standard)) %>%
  select(feature.ID, authentic.standard, authentic.standard.smiles)
if (nrow(level1.annotation) >= 1) {
names(level1.annotation) <- c("feature.ID", "annotation", "smiles")
level1.annotation$confidence.score <- "1"
level1.annotation$confidence.level <- "1"
level1.annotation$annotation.type <- "authentic.standard"
level1.annotation$superclass <- NA
}
not.level1.annotation <- filter(propagation.df, !feature.ID %in% level1.annotation$feature.ID)

gnps <- filter(not.level1.annotation, !is.na(gnps.compound.name)) %>%
  select(feature.ID, gnps.compound.name, gnps.smiles, gnps.cosine.score, gnps.NPC.superclass)
names(gnps) <- c("feature.ID", "annotation", "smiles", "confidence.score", "superclass")
gnps$confidence.level <- "2"
gnps$annotation.type <- "gnps"
not.gnps <- filter(not.level1.annotation, !feature.ID %in% gnps$feature.ID)

csi <- filter(not.gnps, !is.na(csi.compound.name) | csi.compound.name != "null") %>%
  filter(csi.confidence.score > csi.prob) %>%
  select(feature.ID, csi.compound.name, csi.smiles, csi.confidence.score, canopus.NPC.superclass)
names(csi) <- c("feature.ID", "annotation", "smiles", "confidence.score", "superclass")
if (nrow(csi) > 0) {
csi$confidence.level <- "2"
csi$annotation.type <- "csi"
}
not.csi <- filter(not.gnps, !feature.ID %in% csi$feature.ID)

ms2q.real <- filter(not.csi, !is.na(ms2query.analogue.compound.name)) %>%
  filter(ms2query.score > ms2query.prob) %>%
  filter(ms2query.mzdiff <= 0.001) %>%
  select(feature.ID, ms2query.analogue.compound.name, ms2query.smiles, ms2query.score, ms2query.NPC.superclass)
names(ms2q.real) <- c("feature.ID", "annotation", "smiles", "confidence.score", "superclass")
ms2q.real$confidence.level <- "2"
ms2q.real$annotation.type <- "ms2query"
not.real.ms2q <- filter(not.csi, !feature.ID %in% ms2q.real$feature.ID)

ms2q <- filter(not.real.ms2q, !is.na(ms2query.analogue.compound.name)) %>%
  filter(ms2query.score > ms2query.prob) %>%
  select(feature.ID, ms2query.analogue.compound.name, ms2query.smiles, ms2query.score, ms2query.NPC.superclass)
names(ms2q) <- c("feature.ID", "annotation", "smiles", "confidence.score", "superclass")
ms2q$annotation <- paste0("Analogue of ", ms2q$annotation)
ms2q$confidence.level <- "3"
ms2q$annotation.type <- "ms2query analogue"
not.ms2q <- filter(not.real.ms2q, !feature.ID %in% ms2q$feature.ID)

cytoscape.annotations <- rbind(level1.annotation, gnps)
cytoscape.annotations <- rbind(cytoscape.annotations, csi)
cytoscape.annotations <- rbind(cytoscape.annotations, ms2q.real)
cytoscape.annotations <- rbind(cytoscape.annotations, ms2q)

final.annotation.df <-  left_join(propagation.df, cytoscape.annotations, by = "feature.ID")
names(final.annotation.df)[names(final.annotation.df) == "annotation"] <- "Best.Annotation"
names(final.annotation.df)[names(final.annotation.df) == "smiles"] <- "Best.Annotation.Smiles"
names(final.annotation.df)[names(final.annotation.df) == "confidence.score"] <- "Best.Annotation.Confidence.Score"
names(final.annotation.df)[names(final.annotation.df) == "annotation.type"] <- "Best.Annotation.Type"
names(final.annotation.df)[names(final.annotation.df) == "confidence.level"] <- "Best.Annotation.Confidence.Level"
names(final.annotation.df)[names(final.annotation.df) == "superclass"] <- "Best.Annotation.Compound.Class"

names(cytoscape.annotations) <- c("shared.name", "library_compound_name2")

cytoscape.annotations$library_compound_name2 <- ifelse(
  is.na(cytoscape.annotations$library_compound_name2) | cytoscape.annotations$library_compound_name2 == "NA", 
  "", 
  cytoscape.annotations$library_compound_name2
)
names(cytoscape.annotations) <- c("shared.name", "library.compound.name.2", "smiles", "confidence.score", "confidence.level", "annotation.type", "compound.class")
cytoscape.annotations$shared.name <- as.numeric(cytoscape.annotations$shared.name)

cytoscape <- cytoscape %>%
  full_join(cytoscape.annotations, by = "shared.name")
write.csv(cytoscape, paste0(folder, "/cytoscape-v2.csv"))

#Create the samples.df containing feature.usi and peak area.
if ("Samples" %in% names(sample.data2)) {
  sample.data2 <- sample.data2 %>% select(-Samples)
}

# Convert to long format
long_df <- sample.data2 %>%
  pivot_longer(
    cols = -feature.ID,
    names_to = "samples",
    values_to = "area"
  )

samples.df <- final.annotation.df %>%
  select(feature.ID, feature.usi)
samples.df$feature.ID <- as.numeric(samples.df$feature.ID)

samples.df <- long_df %>%
  full_join(samples.df, by = "feature.ID")

#More confidence.level assignments before finalising file.
for (i in 1:nrow(final.annotation.df)) {
  if (is.na(final.annotation.df$Best.Annotation.Confidence.Level[i]) && 
      (!is.na(final.annotation.df$Probable.Analogue.Of[i]))) {
    final.annotation.df$Best.Annotation[i] <- paste0("Probable Analogue of ", final.annotation.df$Probable.Analogue.Of[i])
    final.annotation.df$Best.Annotation.Type[i] <- "Propagated Analogue"
    final.annotation.df$Best.Annotation.Confidence.Level[i] <- "3"
    final.annotation.df$Best.Annotation.Compound.Class[i] <- final.annotation.df$Propagated.Annotation.Class[i]
  }
}

for (i in 1:nrow(final.annotation.df)) {
  if (is.na(final.annotation.df$Best.Annotation.Confidence.Level[i]) && 
      (!is.na(final.annotation.df$ms2query.NPC.pathway[i]) || !is.na(final.annotation.df$canopus.NPC.pathway[i])) &&
      (!is.na(final.annotation.df$ms2query.NPC.pathway[i]) && final.annotation.df$ms2query.NPC.pathway[i] != "None")) { # Check if NOT NA before comparing
    final.annotation.df$Best.Annotation[i] <- paste0("Predicted NPC Pathway: ", final.annotation.df$canopus.NPC.pathway[i])
    final.annotation.df$Best.Annotation.Confidence.Score[i] <- final.annotation.df$canopus.NPC.pathway.probability[i]
    final.annotation.df$Best.Annotation.Confidence.Level[i] <- "4"
    final.annotation.df$Best.Annotation.Type[i] <- "canopus"
    final.annotation.df$Best.Annotation.Compound.Class[i] <- final.annotation.df$canopus.NPC.superclass[i]
  }
}

for (i in 1:nrow(final.annotation.df)) {
  if (is.na(final.annotation.df$Best.Annotation.Confidence.Level[i])) {
    final.annotation.df$Best.Annotation.Confidence.Level[i] <- "5"
  }
}

final.annotation.df <- fix_compound_names(final.annotation.df, "Best.Annotation")  ##df, and column to fix

## 16. Redundancy reduction

dataset <- final.annotation.df %>%
  filter(!is.na(Best.Annotation.Smiles) & Best.Annotation.Smiles != "N/A") # Changed OR to AND
dataset$rt <- as.numeric(dataset$rt)
dataset$redundant <- FALSE

other.data <- final.annotation.df %>%
  filter(is.na(Best.Annotation.Smiles) | Best.Annotation.Smiles == "N/A")

dataset <- redundancy_fixer(dataset, column_to_check = "Best.Annotation.Smiles")

## 17. Writing files to disk      
################################
  ##ANNOTATION FILE      
final.annotation.df2 <- dataset %>% 
  filter(!redundant) %>%
  select(-redundant) %>% 
  rbind(other.data)        
        
write.csv(final.annotation.df2, paste0(folder, "/final-annotation-df.csv"))
write.csv(final.annotation.df2, paste0("Y:/MA_BPA_Microbiome/Dataset-Annotations/", dataset.id, ".csv"))

##################################
  ##SAMPLES DF

    output_file <- paste0(folder, "/", dataset.id, "-samples-df.csv")
    chunk_size <- 50000 # Process 50,000 rows at a time
    total_rows <- nrow(samples.df)
    num_chunks <- ceiling(total_rows / chunk_size)
    
    # Write header row once
    write.table(samples.df[0, ], file = output_file, sep = ",",
                col.names = TRUE, row.names = FALSE, append = FALSE,
                qmethod = "double") # Ensure consistent quoting
    
    cat("Starting to write to CSV...\n")
    
    progress_bar <- txtProgressBar(min = 0, max = num_chunks, style = 3)
    
    for (i in 1:num_chunks) {
      start_row <- (i - 1) * chunk_size + 1
      end_row <- min(i * chunk_size, total_rows)
      
      chunk_data <- samples.df[start_row:end_row, ]
      
      # Write chunk to file, appending to existing content
      write.table(chunk_data, file = output_file, sep = ",",
                  col.names = FALSE, row.names = FALSE, append = TRUE,
                  qmethod = "double")
      
      # Update progress bar
      setTxtProgressBar(progress_bar, i)
      
    }    ##!!PLEASE WAIT!! Samples.df takes a while to write to csv.
    
    close(progress_bar)
    cat(sprintf("\nCSV file saved to: %s\n", output_file))

#Now extracting top 10 features per sample
##First: To Combine Annotations and Samples
final.annotation.df3 <- final.annotation.df2 %>%
  select(feature.ID, mz, Best.Annotation, Best.Annotation.Smiles, Best.Annotation.Confidence.Level)

#Append to Samples
samples.df.with.annotations <- samples.df %>%
  full_join(final.annotation.df3, by = "feature.ID")

##Second: Optional df to extract top 10 best features representative of a sample.
samples.df.with.annotations 

top_10_features <- samples.df.with.annotations %>%
  group_by(samples) %>%
  top_n(10, area) %>%
  ungroup()

write.csv(top_10_features, paste0(folder, "/", dataset.id, "-top-10-features.csv"))

#Script has finished processing