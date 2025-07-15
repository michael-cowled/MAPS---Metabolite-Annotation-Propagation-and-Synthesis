#### ANNOTATION TABLE SCRIPT FOR AUSTRALIAN HUMAN GUT METABOLOME DATABASE ####

#!# Run 1-Functions first! #!!#

#------------------------------------------------------------------------------#
###---START OF USER-FED INFORMATION---###

# Dataset ID: From Data Management Plan:

dataset.id <- "HGMD_0070"     #####Change to dataset of interest#####

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


## 1. check_and_install
# Function to check, install, and load required packages
check_and_install <- function(packages, github_packages = list()) {
  # Install 'remotes' if needed for GitHub installs
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (pkg %in% names(github_packages)) {
        message(paste("Installing", pkg, "from GitHub:", github_packages[[pkg]]))
        remotes::install_github(github_packages[[pkg]])
      } else {
        message(paste("Installing", pkg, "from CRAN"))
        install.packages(pkg, dependencies = TRUE)
      }
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
}

required_packages <- c("MAPS.Package", "dplyr", "tidyr", "stringr", "readr", 
                       "reshape2", "ggplot2", "svglite", "readxl", "data.table", 
                       "openxlsx", "tidyverse", "rvest", "jsonlite", "xml2")

github_packages <- list("MAPS.Package" = "michael-cowled/MAPS-Package-Public")
check_and_install(required_packages, github_packages)

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
folder <- dataset$Processed.Data.Folder[1]

# Run validation and collect paths
paths <- validate_and_get_paths(folder)

# Access each file path as needed
mzmine.data <- paths$mzmine_data
mzmine.annotations <- paths$mzmine_annotations
canopus.data <- paths$canopus_data
csi.data <- paths$csi_data
zodiac.data <- paths$zodiac_data
ms2query.data <- paths$ms2query_data
cytoscape <- paths$cytoscape

## 5. Load in MZMINE data
mzmine.annotations <- read.csv(mzmine.annotations) %>%
  # First, ensure distinct compound_name per id
  distinct(id, compound_name, .keep_all = TRUE) %>%
  group_by(id) %>%
  # Remove all rows with the lowest confidence.score per id
  filter(score > min(score) | n() == 1) %>%
  ungroup() %>%
  # Now, reduce to only one row per compound_name (highest confidence.score overall)
  group_by(compound_name) %>%
  filter(score == max(score)) %>%
  slice(1) %>%  # If there's a tie on confidence.score, pick the first row arbitrarily
  ungroup()

#Calculating ID probability of level 1 annotations
mzmine.annotations$mzmine.id.prob <- NA
#Standardisation of compound names (retrieving from a locally stored cache/pubchem)
# --- Main Processing Loop ---
mzmine.annotations <- standardise_annotation(mzmine.annotations, "compound_name", "smiles")
write.csv(cid_cache_df, cid_cache_file, row.names = FALSE)

# --- Populate from cache ---
for (i in seq_len(nrow(mzmine.annotations))) {
  current_name <- mzmine.annotations$compound_name[i]
  row_match <- cid_cache_df[cid_cache_df$LookupName == current_name, ]
  
  if (nrow(row_match) > 0) {
    mzmine.annotations$CID[i] <- row_match$CID[1]
    mzmine.annotations$smiles[i] <- row_match$SMILES[1]
    mzmine.annotations$mol_formula[i] <- row_match$MolecularFormula[1]
    mzmine.annotations$MonoisotopicMass[i] <- row_match$MonoisotopicMass[1]
  }
}

mzmine.annotations.final <- mzmine.annotations %>%
  group_by(id) %>%
  # Add mzmine.id.prob = 1 / number of rows per id
  mutate(mzmine.id.prob = 1 / n()) %>%
  # Keep only the row(s) with the highest confidence.score per id
  filter(score == max(score)) %>%
  slice(1) %>%  # In case of ties, keep one arbitrarily
  ungroup() %>%
  # Keep only relevant columns:
  select(id, compound_name, score, smiles, mzmine.id.prob, mol_formula, CID, MonoisotopicMass)
names(mzmine.annotations.final) <- c('feature.ID', "compound.name", "confidence.score", 
                                     "smiles", "id.prob", "molecular.formula", "CID", "MonoisotopicMass")
mzmine.annotations.final$feature.ID <- as.numeric(mzmine.annotations.final$feature.ID)
mzmine.annotations.final$confidence.level <- "1"
mzmine.annotations.final$annotation.type <- "authentic standard"

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
names(gnps.annotation.data) <- c("feature.ID", "library.name", "confidence.score", 
                                 "mz.diff.ppm", "gnps.shared.peaks", "compound.name", 
                                 "smiles", "library.quality", "NPC.superclass", "NPC.pathway", "gnps.library.usi")
gnps.annotation.data$gnps.library.usi <- str_replace(
  gnps.annotation.data$gnps.library.usi,
  "mzspec:GNPS:GNPS-LIBRARY:(.*)",
  "mzspec:GNPS:GNPS-LIBRARY:accession:\\1"
)

gnps.annotation.data$gnps.in.silico.bile.acid.info <- NA
gnps.annotation.data$gnps.in.silico.bile.acid.info[grepl("Candidate ", gnps.annotation.data$compound.name)] <- 
  gnps.annotation.data$compound.name[grepl("Candidate ", gnps.annotation.data$compound.name)]
gnps.annotation.data <- gnps.annotation.data[, c(1, 6, 7, 3:5, 2, 8, 10, 9, 11, 12)] #Reorders columns
gnps.annotation.data <- gnps.annotation.data %>%
  filter(mz.diff.ppm <= 5) %>%
  filter(confidence.score >= gnps.prob)

gnps.cluster.data <- read_tsv(paste0("https://gnps2.org/resultfile?task=", gnps.task.id, "&file=nf_output/networking/clustersummary_with_network.tsv"))
gnps.cluster.data <- select(gnps.cluster.data,  'cluster index', 'component')
names(gnps.cluster.data) <- c('feature.ID', "gnps.cluster.ID")
gnps.cluster.pairs <- read_tsv(paste0("https://gnps2.org/resultfile?task=", gnps.task.id, "&file=nf_output/networking/filtered_pairs.tsv"))
gnps.annotation.data$annotation.type <- "gnps"

#Create level 2, and tidy
gnps.data.lv2 <- filter(gnps.annotation.data, library.quality != "Insilico")
gnps.data.lv2$confidence.level <- "2"
gnps.data.lv2 <- fix_compound_names(gnps.data.lv2, "compound.name")

#Create level 3, and tidy
gnps.data.lv3 <- filter(gnps.annotation.data, library.quality == "Insilico")
gnps.data.lv3$confidence.level <- "3"
gnps.data.lv3 <- gnps.data.lv3 %>%
  mutate(
    compound.name = ifelse(
      grepl("Candidate ", compound.name),
      sub("\\s*\\(.*$", "", compound.name),
      compound.name
    )
  )
gnps.data.lv3$gnps.in.silico.bile.acid.info <- gsub('"', "", gnps.data.lv3$gnps.in.silico.bile.acid.info)     
gnps.data.lv3$compound.name <- gsub('"', "", gnps.data.lv3$compound.name)
gnps.data.lv3$CID <- NA
gnps.data.lv3$MolecularFormula <- NA
gnps.data.lv3$MonoisotopicMass <- NA
gnps.data.lv3$id.prob <- NA

## 7. Load in MS2QUERY data
ms2query.data <- read.csv(ms2query.data)
ms2query.data  <- ms2query.data[, c(2, 3, 4, 7, 8, 10, 17, 18)]
names(ms2query.data) <- c("confidence.score", "mz.diff", "precursor_mz", 
                          "compound.name", "smiles", 
                          'feature.ID', "NPC.superclass", "NPC.pathway")
ms2query.data <- fix_compound_names(ms2query.data, "compound.name")  ##df, and column to fix
ms2query.data <- mutate(ms2query.data, mz.diff.ppm = mz.diff/precursor_mz *1000000) %>%
  filter(mz.diff.ppm <= 5)
ms2query.data$annotation.type <- "ms2query"

#Create level 2, and tidy
ms2query.data.lv2 <-ms2query.data %>%
  filter(confidence.score > ms2query.prob) %>%
  filter(mz.diff <= 0.001) %>%
  select(-mz.diff, -precursor_mz)
  
ms2query.data.lv2$confidence.level <- "2"
ms2query.data.lv2$gnps.shared.peaks <- NA
ms2query.data.lv2$library.name <- "ms2query"
ms2query.data.lv2$library.quality <- "ms2query"
ms2query.data.lv2$gnps.library.usi <- NA
ms2query.data.lv2$gnps.in.silico.bile.acid.info <- NA

#Create level 3, and tidy
ms2query.data.lv3 <- ms2query.data %>%
  filter(confidence.score > ms2query.prob) %>%
  filter(mz.diff >= 0.001) %>%
  filter(!feature.ID %in% mzmine.annotations.final$feature.ID) #removes any annotations that conflict with level 1 annotations
ms2query.data.lv3$confidence.level <- "3"

## 8. Standardise level 2 annotations and compute ID probability
lv2.annotations <- rbind(gnps.data.lv2, ms2query.data.lv2)  ##bind the two sets of annotations together

lv2.annotations$CID <- NA
lv2.annotations$MolecularFormula<- NA
lv2.annotations$MonoisotopicMass <- NA
lv2.annotations <- deduplicate_data(lv2.annotations, compound.name, confidence.score) #Pre-standardisation filtering for duplicates
lv2.annotations<- standardise_annotation(lv2.annotations, "compound.name", "smiles") #Standardisation of naming using pubchem, remove for time-savings
lv2.annotations <- deduplicate_data(lv2.annotations, compound.name, confidence.score) #Post-standardisation filtering of duplicates
write.csv(cid_cache_df, cid_cache_file, row.names = FALSE)

#Removal of lv1 annotations
lv2.annotations <- filter_data_by_standards(lv2.annotations, mzmine.data, "compound.name", "feature.ID", rt.tol)

#Computing of ID prob.
lv2.annotations <- lv2.annotations %>%
  group_by(feature.ID) %>%
  mutate(
    n_above_thresh = sum(confidence.score >= ms2query.prob),
    id.prob = ifelse(confidence.score >= ms2query.prob & n_above_thresh > 0, 1 / n_above_thresh, 0)
  ) %>%
  ungroup() %>%
  select(-n_above_thresh) %>%
  group_by(feature.ID) %>%
  arrange(desc(confidence.score)) %>%
  slice(1) %>%
  ungroup() %>%
  select(-rt)

## 9. Appending lv3 In silico matches from GNPS
lv2.and.lv3.annotations <- lv2.annotations %>%
  rbind(gnps.data.lv3)

## 10. Load in SIRIUS data:      Compatible with v6.1.0 onwards
canopus.data <- read_tsv(canopus.data)
canopus.data <- canopus.data[, c(5:8, 27)]
names(canopus.data) <- c("canopus.NPC.pathway", "canopus.NPC.pathway.probability", 
                         "canopus.NPC.superclass", "canopus.NPC.superclass.probability", 
                         'feature.ID')
canopus.data <- canopus.data %>%
  group_by(feature.ID) %>%
  filter(!(all(canopus.NPC.pathway.probability == 0))) %>%  # Remove groups where all confidence.scores are 0
  filter(canopus.NPC.pathway.probability == max(canopus.NPC.pathway.probability, na.rm = TRUE)) %>%  # Keep the highest confidence.score
  slice(1) %>%  # Arbitrarily keep one if multiple rows have the same non-zero confidence.score
  ungroup()

zodiac.data <- read_tsv(zodiac.data)
zodiac.data  <- zodiac.data[, c(2, 5, 21)]
names(zodiac.data) <- c("zodiac.formula", "zodiac.confidence.score", 'feature.ID')
# Filter and keep only unique feature.ID (in the case where multiple annotations are provided)
zodiac.data <- zodiac.data %>%
  group_by(feature.ID) %>%
  filter(!(all(zodiac.confidence.score == 0))) %>%  # Remove groups where all confidence.scores are 0
  filter(zodiac.confidence.score == max(zodiac.confidence.score, na.rm = TRUE)) %>%  # Keep the highest confidence.score
  slice(1) %>%  # Arbitrarily keep one if multiple rows have the same non-zero confidence.score
  ungroup()

csi.data <- read_tsv(csi.data)
csi.data <- csi.data[, c(3, 14, 15, 25)]
names(csi.data) <- c("confidence.score", "compound.name", "smiles", 'feature.ID')
csi.data <- csi.data[, c(2, 1, 3:ncol(csi.data))] #Swap cols 1 and 2
# Initial fixes
csi.data$confidence.score[csi.data$confidence.score == -Inf] <- 0
csi.data$confidence.score <- as.numeric(csi.data$confidence.score)
csi.data <- filter(csi.data, confidence.score >= csi.prob)
csi.data$compound.name[grepl("Solaparnaine", csi.data$compound.name, ignore.case = TRUE)] <- "Solaparnaine" ##troublesome case
csi.data$compound.name[grepl("Spectalinine", csi.data$compound.name, ignore.case = TRUE)] <- "(-)-Spectalinine" ##troublesome case

## 11. Appending lv3 In silico matches from CSI:FingerID
# Calculation of Identification Probability (Wishart paper)
csi.data <- csi.data %>%
  mutate(id.prob = 0)
csi.data <- csi.data %>%
  group_by(feature.ID) %>%
  mutate(
    n_above_thresh = sum(confidence.score >= csi.prob),
    id.prob = ifelse(confidence.score >= csi.prob & n_above_thresh > 0, 1 / n_above_thresh, 0)
  ) %>%
  ungroup() %>%
  select(-n_above_thresh) %>%
  group_by(feature.ID) %>%
  arrange(desc(confidence.score)) %>%
  slice(1) %>%
  ungroup()

#Filter for unique entries to csi.data not already in lv2 and lv3 annotations
unique_in_csi <- setdiff(csi.data$feature.ID, lv2.and.lv3.annotations$feature.ID)
csi.data <- filter(csi.data, feature.ID %in% unique_in_csi)

#################Evenutally add in standardisation and removal of lv1 AND lv2 from here ---------------------------------------------------------------------------------------## TO DO
#[in similar style to what was done for gnps. Do at this point after reducing, and prob. is calculated diff. i.e. estimate #of hits, and remove by 1.]

# Identify columns present in lv2.and.lv3.annotations but not in csi.data
missing_cols <- colnames(lv2.and.lv3.annotations)[!colnames(lv2.and.lv3.annotations) %in% colnames(csi.data)]
for (col in missing_cols) { # Add the missing columns to csi.data and fill with NA
  csi.data[[col]] <- NA
}
csi.data$annotation.type <- "CSI:FingerID"
csi.data$confidence.level <- "3"

lv2.and.lv3.annotations <- lv2.and.lv3.annotations %>%
  rbind(csi.data)

## 12. Appending lv3 analogues from ms2query
ms2query.data.lv3$compound.name <- paste0("Analogue of ", ms2query.data.lv3$compound.name)
#Filter for unique entries to csi.data not already in lv2 and lv3 annotations
unique_in_ms2query <- setdiff(ms2query.data.lv3$feature.ID, lv2.and.lv3.annotations$feature.ID)
ms2query.data.lv3 <- filter(ms2query.data.lv3, feature.ID %in% unique_in_ms2query) %>%
  select(-mz.diff, -precursor_mz)

# Identify columns present in lv2.and.lv3.annotations but not in csi.data
missing_cols <- colnames(lv2.and.lv3.annotations)[!colnames(lv2.and.lv3.annotations) %in% colnames(ms2query.data.lv3)]
for (col in missing_cols) { # Add the missing columns to ms2query.data.lv3 and fill with NA
  ms2query.data.lv3[[col]] <- NA
}

lv2.and.lv3.annotations <- lv2.and.lv3.annotations %>%
  rbind(ms2query.data.lv3)

## 13. Appending lv1 annotations and all other annotations
missing_cols <- colnames(lv2.and.lv3.annotations)[!colnames(lv2.and.lv3.annotations) %in% colnames(mzmine.annotations.final)]
for (col in missing_cols) { # Add the missing columns to ms2query.data.lv3 and fill with NA
  mzmine.annotations.final[[col]] <- NA
}
mzmine.annotations.final$confidence.level <- "1"
mzmine.annotations.final <- select(mzmine.annotations.final, -molecular.formula)

lv1.lv2.lv3.annotations <- lv2.and.lv3.annotations %>%
  rbind(mzmine.annotations.final)

#Merge annotations into one big table                                     
full.annotation.data <- mzmine.data[,1:4] %>%
  full_join(lv2.and.lv3.annotations, by = "feature.ID") %>%
  full_join(canopus.data, by = "feature.ID") %>%
  full_join(zodiac.data, by = "feature.ID") %>%
  full_join(gnps.cluster.data, by = "feature.ID")

full.annotation.data$feature.usi <- paste0("mzspec:GNPS2:TASK-", gnps.task.id, "-nf_output/clustering/spectra_reformatted.mgf:scan:", full.annotation.data$feature.ID) #Adds unique spectra identifiers

## 14. Propagation of annotations
# Initialize columns in propagated.annotation.data
propagated.annotation.data <- full.annotation.data %>%
  mutate(
    Probable.Analogue.Of = NA,
    Propagated.Feature.ID = NA,
    Propagated.Annotation.Type = NA, # Initialize the new column
    Propagated.Annotation.Class = NA # This is already here
  )

# Identify the unknown features
na.rows <- filter(propagated.annotation.data, is.na(compound.name))
na.feature.ids <- na.rows$feature.ID

# Iterate over na.feature.ids to populate propagated.annotation.data
for (i in na.feature.ids) {
  paired_values <- paired_feature_finder(i)
  print(paste("Paired values for feature ID", i, ":", paste(paired_values, collapse = ", ")))
  
  result_data <- list(value = NA, column = NA, superclass = NA)
  selected_paired_value <- NA
  
  for (value in paired_values) {
    result_data <- get_result(value, full.annotation.data)
    if (!is.na(result_data$value)) {
      selected_paired_value <- value
      break
    }
  }
  
  # If you want to ensure annotation.type specifically:
  if (!is.na(selected_paired_value)) {
    annotation_type_value <- full.annotation.data %>%
      filter(feature.ID == selected_paired_value) %>%
      pull(annotation.type) %>%
      first()
    
    result_data$column <- annotation_type_value
  }
  
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
      Propagated.Annotation.Class = case_when(
        feature.ID == i ~ result_data$superclass,
        TRUE ~ Propagated.Annotation.Class
      )
    )
}

# Calculating the mz difference of the analogue to the propagated feature:
# Convert mz column to numeric if it's not already
if (!is.numeric(full.annotation.data$mz)) {
  full.annotation.data$mz <- as.numeric(full.annotation.data$mz)
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
    mz_value <- full.annotation.data$mz[full.annotation.data$feature.ID == feature_id]
    propagated_mz_value <- full.annotation.data$mz[full.annotation.data$feature.ID == propagated_feature_id]
    
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

##Append level 4 and 5 annotations
propagated.annotation.data <- propagated.annotation.data %>%
  mutate(
    confidence.level = as.character(confidence.level),  # ensure editable
    level4_mask = is.na(compound.name) & !is.na(canopus.NPC.pathway),  # condition mask
    compound.name = ifelse(
      level4_mask,
      paste0("Predicted NPC Pathway: ", canopus.NPC.pathway),
      compound.name
    ),
    confidence.level = ifelse(
      level4_mask,
      "4",
      confidence.level
    ),
    NPC.pathway = canopus.NPC.pathway,
    NPC.superclass = canopus.NPC.superclass,
    annotation.type = ifelse(
      level4_mask,
      "canopus",
      annotation.type
    )
  ) %>%
  select(-level4_mask) %>%
  mutate(
    confidence.level = ifelse(is.na(confidence.level), "5", confidence.level)
  )

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

## 14. Cytoscape Appending
#Simultaneously creates a tidied cytoscape file with the best annotation (including analogues)
#Represented as the output.

cytoscape <- read.csv(paste0(folder, "/gnps/cytoscape.csv"))

cytoscape.annotations <- select(propagated.annotation.data.with.samples, feature.ID, compound.name)
names(cytoscape.annotations) <- c("shared.name", "library_compound_name2")

cytoscape.annotations$library_compound_name2 <- ifelse(
  is.na(cytoscape.annotations$library_compound_name2) | cytoscape.annotations$library_compound_name2 == "NA", 
  "", 
  cytoscape.annotations$library_compound_name2
)
cytoscape.annotations$shared.name <- as.numeric(cytoscape.annotations$shared.name)
cytoscape <- cytoscape %>%
  full_join(cytoscape.annotations, by = "shared.name")
colnames(cytoscape)[colnames(cytoscape) == 'shared.name'] <- 'shared name'

write.csv(cytoscape, paste0(folder, "/cytoscape-v2.csv"), row.names = FALSE)

## 15. Collapsing Ion Identity Networks
#performed in such a way to retain the best annotation (if multiple)

#A new editing df where features are sequentially removed:
final.annotation.df <- propagated.annotation.data.with.samples

#A new df with "good" annotations extracted in a similar manner to unpropagated version.
new.summary.annotation.data <- propagated.annotation.data.with.samples

iin.features <- filter(new.summary.annotation.data, !is.na(ion.identity.ID))

results <- process_all_features(data, "ion.identity.ID", c("feature1", "feature2"))

# Print results for each ID
for (id in names(results)) {
  cat("Results for ion.identity.ID =", id, ":\n")
  print(results[[id]])
}

# Update final.annotation.df by removing rows based on results and retaining NA IDs
final.annotation.df <- update_data_frame(final.annotation.df, results, "ion.identity.ID") %>%
  select(-ion.identity.ID, -Probable.Analogue.Of)


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
  select(feature.ID, feature.usi, compound.name, smiles)
samples.df$feature.ID <- as.numeric(samples.df$feature.ID)

samples.df <- long_df %>%
  full_join(samples.df, by = "feature.ID")

final.annotation.df <- fix_compound_names(final.annotation.df, "compound.name")  ##df, and column to fix

## 16. Redundancy reduction

dataset <- final.annotation.df %>%
  filter(!is.na(smiles) & smiles != "N/A") # Changed OR to AND
dataset$rt <- as.numeric(dataset$rt)
dataset$redundant <- FALSE

other.data <- final.annotation.df %>%
  filter(is.na(smiles) | smiles == "N/A")

dataset <- redundancy_fixer(dataset, column_to_check = "smiles")

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
  select(feature.ID, mz, compound.name, smiles, confidence.level)

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