#### ANNOTATION TABLE SCRIPT FOR AUSTRALIAN HUMAN GUT METABOLOME DATABASE ####

###---START OF USER-FED INFORMATION---###

# Dataset ID: From Data Management Plan:
dataset.id <- "HGMD_0050"     #####Change to dataset of interest#####

# Specify the path to the Data Management Plan (unless moved to Mediaflux then should be consistent)
excel_file <- "C:/Users/mcowled/The University of Melbourne/Gut microbiome group - Documents/Data Management and Analysis/HGM - Data Management System.xlsx"

# Annotation Acceptance Probabilities (Defaults for Exploratory Analyses)
# Increase stringency if the application demands it.

canopus.prob <- 0.7           # Default is 0.7
csi.prob <- 0.8863            # Default is 0.8863
ms2query.prob <- 0.63         # Default is 0.63

###---END OF USER-FED INFORMATION---###

# Loading libraries

library(ggplot2)
library(reshape2) 
library(dplyr)
library(tidyr)
library(readr)
library(svglite)
library(readxl)

##-Updating-Metadata-##

sheet_names <- excel_sheets(excel_file) # Get the sheet names
for (sheet in sheet_names) {
  data <- read_excel(excel_file, sheet = sheet)
  write.csv(data, file = paste0(sheet, ".csv"), row.names = FALSE)
}

#GNPS2 Task ID: See URL for job: https://gnps2.org/status?task=TASKID 
dataset <- read.csv("D - Dataset.csv") %>%
  filter(HGMD.ID == dataset.id)
gnps.task.id <- dataset$gnps.task.ID[1] #TASKID = unique set of random numbers and letters

##Import Libraries

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
required_packages <- c("dplyr", "tidyr", "stringr", "readr") # Check, install, and load required packages
check_and_install(required_packages) # Now your packages are installed and loaded

## Import annotation tables to merge - !!! Check spelling !!!

# mzmine - works for v4.4.0 and specific HGMA batch file
mzmine.data <- read.csv("ms1-and-ms2.csv") # Derived from Export to CSV file (modular)
sample.data <- mzmine.data # A copy to be used for different processing later
mzmine.data <- select(mzmine.data, "id", "rt", "mz", "ion_identities.iin_id", "spectral_db_matches.compound_name") %>%
  mutate(across(everything(), as.character)) %>%
  mutate(across(everything(), ~ na_if(., "")))
names(mzmine.data) <- c('feature.ID', "rt", "mz", "ion.identity.ID", "authentic.standard")
mzmine.data$feature.ID <- as.numeric(mzmine.data$feature.ID)
mzmine.data$authentic.standard.smiles <- NA                         #Set up to scrape SMILES later

# GNPS2 - works for v0.1.2, no metadata required
gnps.annotation.data <- read_tsv(paste0("https://gnps2.org/resultfile?task=", gnps.task.id, "&file=nf_output/library/merged_results_with_gnps.tsv"))
gnps.annotation.data <- gnps.annotation.data[, c(2, 4, 5, 8, 9, 15, 27, 35, 43, 45, 46)]
names(gnps.annotation.data) <- c("feature.ID", "gnps.library.name", "gnps.cosine.score", 
                                 "gnps.diff.ppm", "gnps.shared.peaks", "gnps.compound.name", 
                                 "gnps.smiles", "gnps.library.quality", "gnps.NPC.superclass", "gnps.NPC.pathway", "library.usi")
gnps.annotation.data <- gnps.annotation.data[, c(1, 6, 7, 3:5, 2, 8, 10, 9, 11)] #Reorders columns
gnps.cluster.data <- read_tsv(paste0("https://gnps2.org/resultfile?task=", gnps.task.id, "&file=nf_output/networking/clustersummary_with_network.tsv"))
gnps.cluster.data <- select(gnps.cluster.data,  'cluster index', 'component')
names(gnps.cluster.data) <- c('feature.ID', "gnps.cluster.ID")
gnps.data <- gnps.cluster.data %>%
  full_join(gnps.annotation.data, by = "feature.ID")
gnps.data$feature.usi <- paste0("mzspec:GNPS2:TASK-", gnps.task.id, "-nf_output/clustering/spectra_reformatted.mgf:scan:", gnps.data$feature.ID) #Adds unique spectra identifiers
gnps.cluster.pairs <- read_tsv(paste0("https://gnps2.org/resultfile?task=", gnps.task.id, "&file=nf_output/networking/filtered_pairs.tsv"))

        ###Tidying GNPS naming conventions
        # Create a new column 'gnps.in.silico.bile.acid.info'
        gnps.data$gnps.in.silico.bile.acid.info <- NA
        
        # Copy values containing "Candidate " to the new column
        gnps.data$gnps.in.silico.bile.acid.info[grepl("Candidate ", gnps.data$gnps.compound.name)] <- 
          gnps.data$gnps.compound.name[grepl("Candidate ", gnps.data$gnps.compound.name)]
        
        fix_compound_names <- function(df) {
          # Convert to data.table for faster operations
          df <- data.table::as.data.table(df)
          
          # Define patterns and replacements
          patterns <- c(
            "\\.(alpha|beta|gamma|omega)\\.-", "\\.(DELTA)\\.", "[[:space:]]*\\[IIN-based(?:[[:space:]]+on)?[[:space:]]*:[[:space:]]*[^)]*\\]", 
            "NCG[[:alnum:]-]+_[[:alnum:]]+_", "NCG[^\\s!]+[\\s!]?", 
            "Spectral Match to ", " from NIST14",  
            "40.0 eV", "70.0 eV", "50.0 eV", ".00 eV", '"', "CollisionEnergy",  "\\|.*", 
            "_", "\\[", "\\]", "\\{", "\\}", " :102040", " :205060", " - ",
            "30.0 eV", "Massbank:[[:alnum:]]+\\s? ", "MassbankEU:[[:alnum:]]+\\s? ", 
            "MoNA:[[:alnum:]]+\\s? ", " (Chimeric precursor selection)", "; \\(M\\+.*" 
            
            
          )
          replacements <- c(
            "\\1-", "delta-", "", "", "", "", "", "", "", "", "", "",  "", "", 
            " ", "(", ")", "(", ")", "", "",  "", "", "", "", "", "", "", "" 
          )
          
          # Columns to apply the fixes
          columns <- "gnps.compound.name"
          
          # Loop through columns and apply fixes
          for (col in columns) {
            for (i in seq_along(patterns)) {
              df[grepl(patterns[i], get(col), ignore.case = TRUE), (col) := gsub(patterns[i], replacements[i], get(col))]
            }
          }
          
          return(as.data.frame(df))
        }
        
        # Apply the function to your data frame
        gnps.data <- fix_compound_names(gnps.data)
        
        gnps.data <- gnps.data %>%
          mutate(
            gnps.compound.name = ifelse(
              grepl("Candidate ", gnps.compound.name),
              sub("\\s*\\(.*$", "", gnps.compound.name),
              gnps.compound.name
            )
          )
        gnps.data$gnps.in.silico.bile.acid.info <- gsub('"', "", gnps.data$gnps.in.silico.bile.acid.info)
        
# SIRIUS:      Compatible with v6.1.0 onwards
canopus.data <- read_tsv('canopus_structure_summary.tsv')
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

csi.data <- read_tsv('structure_identifications.tsv')
csi.data <- csi.data[, c(3, 14, 15, 25)]
names(csi.data) <- c("csi.confidence.score", "csi.compound.name", "csi.smiles", 'feature.ID')
csi.data <- csi.data[, c(2, 1, 3:ncol(csi.data))] #Swap cols 1 and 2
# Filter and keep only unique feature.ID (in the case where multiple annotations are provided)
csi.data <- csi.data %>%
  group_by(feature.ID) %>%
  filter(!(all(csi.confidence.score == 0))) %>%  # Remove groups where all scores are 0
  filter(csi.confidence.score == max(csi.confidence.score, na.rm = TRUE)) %>%  # Keep the highest score
  slice(1) %>%  # Arbitrarily keep one if multiple rows have the same non-zero score
  ungroup()
# Convert the `csi.confidence.score` column to numeric, handling errors
csi.data$csi.confidence.score <- as.numeric(csi.data$csi.confidence.score)
# Replace '-Infinity' with 0 in the csi.confidence.score column
csi.data$csi.confidence.score[csi.data$csi.confidence.score == -Inf] <- 0


zodiac.data <- read_tsv('formula_identifications.tsv')
zodiac.data  <- zodiac.data[, c(2, 5, 21)]
names(zodiac.data) <- c("zodiac.formula", "zodiac.score", 'feature.ID')
# Filter and keep only unique feature.ID (in the case where multiple annotations are provided)
zodiac.data <- zodiac.data %>%
  group_by(feature.ID) %>%
  filter(!(all(zodiac.score == 0))) %>%  # Remove groups where all scores are 0
  filter(zodiac.score == max(zodiac.score, na.rm = TRUE)) %>%  # Keep the highest score
  slice(1) %>%  # Arbitrarily keep one if multiple rows have the same non-zero score
  ungroup()

# MS2QUERY:write
ms2query.data <- read.csv("ms2query.csv")
ms2query.data  <- ms2query.data[, c(2, 3, 7, 8, 10, 17, 18)]
names(ms2query.data) <- c("ms2query.score", "ms2query.mzdiff", 
                          "ms2query.analogue.compound.name", "ms2query.smiles", 
                          'feature.ID', "ms2query.NPC.superclass", "ms2query.NPC.pathway")
ms2query.data <- ms2query.data[, c(3, 1, 2, 4, 5, 7, 6)] #Reorders columns

## Merge annotations into one big table                                     
full.annotation.data <- mzmine.data %>%
  full_join(gnps.data, by = "feature.ID") %>%
  full_join(canopus.data, by = "feature.ID") %>%
  full_join(csi.data, by = "feature.ID") %>%
  full_join(zodiac.data, by = "feature.ID") %>%
  full_join(ms2query.data, by = "feature.ID")

#PubChem Scraper - Finds SMILES if missing

library(httr)

# Function to fetch Synonyms from PubChem, handling errors
get_pubchem <- function(name, input, output) { # Typical inputs = 'name', 'cid'; outputs = 'property/SMILES', 'synonyms', 'cids'
  if (is.na(name)) {
    return(NA) 
  }
  
  url <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/", input, "/", 
                URLencode(name), "/", output, "/TXT")
  
  tryCatch({
    response <- GET(url)
    if (http_error(response)) {
      warning(paste("Request failed for", name, "with status code:", status_code(response)))
      return(NA)
    } else {
      Synonyms <- as.character(response)
      return(Synonyms)
    }
  }, error = function(e) {
    warning(paste("Error fetching Synonyms for", name, ":", e))
    return(NA)
  })
}

# Get SMILES based on COMPOUND NAME via CIDS - Authentic Standard

smiles_list <- character(nrow(full.annotation.data))  # Initialize with correct length

for (i in seq_len(nrow(full.annotation.data))) {
  if (!is.na(full.annotation.data$authentic.standard[i]) && 
      (is.na(full.annotation.data$authentic.standard.smiles[i]) || full.annotation.data$authentic.standard.smiles[i] == "N/A")) {
    cid <- get_pubchem(full.annotation.data$authentic.standard[i], "name", "cids")
    cid <- as.character(gsub("\n.*", "", cid))
    print(full.annotation.data$authentic.standard[i])
    print(cid)
    Sys.sleep(0.2)
    
    smiles <- as.character(get_pubchem(cid, "cid", "property/SMILES"))
    smiles <- as.character(gsub("\n.*", "", smiles))
    
    print(smiles)
    
    smiles_list[i] <- smiles  # Store the SMILES in the list
    print(smiles_list[i])
    
    Sys.sleep(0.2)
    print(paste("Processed row:", i))
  } else {
    smiles_list[i] <- full.annotation.data$authentic.standard.smiles[i]  # Keep original value
  }
}

full.annotation.data$authentic.standard.smiles <- smiles_list

# Get SMILES based on COMPOUND NAME via CIDS - GNPS

smiles_list <- character(nrow(full.annotation.data))  # Initialize with correct length

for (i in seq_len(nrow(full.annotation.data))) {
  if (!is.na(full.annotation.data$gnps.compound.name[i]) && 
      (is.na(full.annotation.data$gnps.smiles[i]) || full.annotation.data$gnps.smiles[i] == "N/A") && 
      (!grepl("Candidate", full.annotation.data$gnps.compound.name[i]))) { 
    
    cid <- get_pubchem(full.annotation.data$gnps.compound.name[i], "name", "cids")
    cid <- as.character(gsub("\n.*", "", cid))
    print(full.annotation.data$gnps.compound.name[i])
    print(cid)
    Sys.sleep(0.2)
    
    smiles <- as.character(get_pubchem(cid, "cid", "property/SMILES"))
    smiles <- as.character(gsub("\n.*", "", smiles))
    
    print(smiles)
    
    smiles_list[i] <- smiles  # Store the SMILES in the list
    print(smiles_list[i])
    
    Sys.sleep(0.2)
    print(paste("Processed row:", i))
  } else {
    smiles_list[i] <- full.annotation.data$gnps.smiles[i]  # Keep original value
  }
}

full.annotation.data$gnps.smiles <- smiles_list

###Propagation Table
##Using full.annotation.data, perform propagations

##Create a summary.annotation.data table, that is filtered for all annotations satisfying a certain degree of confidence
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
    Propagated.Annotation.Type = NA  # Initialize the new column
  )

# Function to get the first non-NA value from specified columns and the column name
get_result <- function(paired_value, summary_data) {
  # Define the columns to check for non-NA values
  columns_to_check <- c("gnps.compound.name", "csi.compound.name", 
                        "ms2query.analogue.compound.name")
  
  # Filter the relevant row and select the columns
  data_subset <- summary_data %>%
    filter(feature.ID == paired_value) %>%
    select(all_of(columns_to_check), csi.smiles)
  
  # Find the first non-NA and non-'null' value
  result <- NA
  column_name <- NA
  for (col in columns_to_check) {
    value <- data_subset[[col]]
    if (!is.na(value) && value != 'null') {
      result <- value
      column_name <- col
      break
    }
  }
  
  # If no valid result is found, fallback to 'csi.smiles'
  if (is.na(result) || result == 'null') {
    csi_value <- data_subset$csi.smiles
    if (!is.na(csi_value) && csi_value != 'null') {
      result <- csi_value
      column_name <- "csi.smiles"
    }
  }
  
  # Return the result and corresponding column name
  list(value = result, column = column_name)
}

# Function: Finds the features linked by an edge in a GNPS cluster
paired.feature.finder <- function(ID) {
  filtered_pairs <- gnps.cluster.pairs %>%
    filter(CLUSTERID1 == ID | CLUSTERID2 == ID) %>%
    arrange(desc(Cosine))  # Sort by decreasing Cosine
  
  # Extract paired values
  paired_values <- filtered_pairs %>%
    mutate(paired_value = ifelse(CLUSTERID1 == ID, CLUSTERID2, CLUSTERID1)) %>%
    pull(paired_value)
  
  return(paired_values)
}

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
  paired_values <- paired.feature.finder(i)
  print(paste("Paired values for feature ID", i, ":", paste(paired_values, collapse = ", ")))
  
  # Initialize variables to store results and corresponding paired value
  result_data <- list(value = NA, column = NA)
  selected_paired_value <- NA
  
  for (value in paired_values) {
    result_data <- get_result(value, summary.annotation.data)
    if (!is.na(result_data$value)) {
      selected_paired_value <- value  # Store the paired value used for the result
      break  # Exit the loop if a valid result is found
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
      )
    )
}

#Calculating the mz difference of the analogue to the propagated feature:
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

##Now to append sample.list corresponding to each feature

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

# Replace all positive areas with 1 i.e. presence.absence
colnames(sample.data) <- colnames_sample
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

##Collapsing Ion Identity Networks - performed in such a way to retain the best annotation (if multiple)

#A new editing df where features are sequentially removed:
propagation.df <- propagated.annotation.data.with.samples

#A new df with "good" annotations extracted in a similar manner to unpropagated version.
new.summary.annotation.data <- propagated.annotation.data.with.samples
#Apply the conditions
new.summary.annotation.data$csi.compound.name[new.summary.annotation.data$csi.confidence.score < csi.prob] <- NA
new.summary.annotation.data$csi.smiles[new.summary.annotation.data$csi.confidence.score < csi.prob] <- NA
new.summary.annotation.data$ms2query.analogue.compound.name[new.summary.annotation.data$ms2query.score < ms2query.prob] <- NA

iin.features <- filter(new.summary.annotation.data, !is.na(ion.identity.ID))

# Define a function to process each IIN group
process_iin_group <- function(group_data, feature_columns) {
  # Initialize a vector to track which rows to retain
  rows_to_retain <- logical(nrow(group_data))
  
  for (col in feature_columns) {
    if (any(!is.na(group_data[[col]]))) {
      values <- group_data$feature.ID[!is.na(group_data[[col]])]
      if (length(values) > 1) {
        # Multiple values meet the condition, retain the lowest values
        highest_value <- max(values, na.rm = TRUE)
        rows_to_retain <- group_data$feature.ID %in% values & group_data$feature.ID != highest_value
      } else {
        # Only one value meets the condition, retain this row
        rows_to_retain <- group_data$feature.ID == values
      }
      break  # Exit loop once a condition is met
    }
  }
  
  if (all(!rows_to_retain)) {
    # No rows were marked to retain, so retain the one with the minimum feature.ID for fallback
    rows_to_retain <- group_data$feature.ID == min(group_data$feature.ID, na.rm = TRUE)
  }
  
  return(group_data[rows_to_retain, ])
}

# Define a function to update the main data frame
update_data_frame <- function(df, results, id_column) {
  df_non_na <- df %>% filter(!is.na(!!sym(id_column)))
  results_df <- bind_rows(results)
  
  updated_df <- df_non_na %>% filter(feature.ID %in% results_df$feature.ID)
  df_na <- df %>% filter(is.na(!!sym(id_column)))
  
  final_df <- bind_rows(updated_df, df_na)
  return(final_df)
}

# Main processing function
process_all_features <- function(data, id_column, feature_columns) {
  unique_ids <- unique(data[[id_column]][!is.na(data[[id_column]])])
  results <- lapply(unique_ids, function(id) {
    group_data <- filter(data, !!sym(id_column) == id)
    process_iin_group(group_data, feature_columns)
  })
  names(results) <- as.character(unique_ids)
  return(results)
}

# Print results for each ID
for (id in names(results)) {
  cat("Results for ion.identity.ID =", id, ":\n")
  print(results[[id]])
}

# Update propagation.df by removing rows based on results and retaining NA IDs
propagation.df <- update_data_frame(propagation.df, results, "ion.identity.ID")

####Cytoscape Tidier####

##Uses output from Annotation-Table.R##

cytoscape <- read.csv("cytoscape.csv")

level1.annotation <- filter(propagation.df, !is.na(authentic.standard)) %>%
  select(feature.ID, authentic.standard, authentic.standard.smiles)
names(level1.annotation) <- c("feature.ID", "annotation", "smiles")
level1.annotation$confidence.level <- "1"
not.level1.annotation <- filter(propagation.df, !feature.ID %in% level1.annotation$feature.ID)

gnps <- filter(not.level1.annotation, !is.na(gnps.compound.name)) %>%
  select(feature.ID, gnps.compound.name, gnps.smiles)
names(gnps) <- c("feature.ID", "annotation", "smiles")
gnps$confidence.level <- "2"
not.gnps <- filter(not.level1.annotation, !feature.ID %in% gnps$feature.ID)

csi <- filter(not.gnps, !is.na(csi.compound.name) | csi.compound.name != "null") %>%
  filter(csi.confidence.score > csi.prob) %>%
  select(feature.ID, csi.compound.name, csi.smiles)
names(csi) <- c("feature.ID", "annotation", "smiles")
csi$confidence.level <- "2"
not.csi <- filter(not.gnps, !feature.ID %in% csi$feature.ID)

ms2q.real <- filter(not.csi, !is.na(ms2query.analogue.compound.name)) %>%
  filter(ms2query.score > ms2query.prob) %>%
  filter(ms2query.mzdiff <= 0.001) %>%
  select(feature.ID, ms2query.analogue.compound.name, ms2query.smiles)
names(ms2q.real) <- c("feature.ID", "annotation", "smiles")
ms2q.real$confidence.level <- "2"
not.real.ms2q <- filter(not.csi, !feature.ID %in% ms2q.real$feature.ID)

ms2q <- filter(not.real.ms2q, !is.na(ms2query.analogue.compound.name)) %>%
  filter(ms2query.score > ms2query.prob) %>%
  select(feature.ID, ms2query.analogue.compound.name, ms2query.smiles)
names(ms2q) <- c("feature.ID", "annotation", "smiles")
ms2q$annotation <- paste0("Analogue of ", ms2q$annotation)
ms2q$confidence.level <- "3"
not.ms2q <- filter(not.real.ms2q, !feature.ID %in% ms2q$feature.ID)

cytoscape.annotations <- rbind(level1.annotation, gnps)
cytoscape.annotations <- rbind(cytoscape.annotations, csi)
cytoscape.annotations <- rbind(cytoscape.annotations, ms2q.real)

final.annotation.df <-  left_join(propagation.df, cytoscape.annotations, by = "feature.ID") %>%
  select(-confidence.level)
names(final.annotation.df)[names(final.annotation.df) == "annotation"] <- "Best.Annotation"
names(final.annotation.df)[names(final.annotation.df) == "smiles"] <- "Best.Annotation.Smiles"

cytoscape.annotations <- rbind(cytoscape.annotations, ms2q)

final.annotation.df <-  left_join(final.annotation.df, cytoscape.annotations, by = "feature.ID")
names(final.annotation.df)[names(final.annotation.df) == "annotation"] <- "Best.Annotation.with.Analogues"
names(final.annotation.df)[names(final.annotation.df) == "smiles"] <- "Best.Annotation.with.Analogues.Smiles"

names(cytoscape.annotations) <- c("shared.name", "library_compound_name2")

cytoscape.annotations$library_compound_name2 <- ifelse(
  is.na(cytoscape.annotations$library_compound_name2) | cytoscape.annotations$library_compound_name2 == "NA", 
  "", 
  cytoscape.annotations$library_compound_name2
)
names(cytoscape.annotations) <- c("shared.name", "library.compound.name.2", "smiles", "confidence.level")

cytoscape <- cytoscape %>%
  full_join(cytoscape.annotations, by = "shared.name")
write.csv(cytoscape, "cytoscape-v2.csv")

samples.df <- final.annotation.df %>%
  select(feature.ID, feature.usi, Samples) %>%
  separate_rows(Samples, sep = "; ")
write.csv(samples.df, "samples-df.csv")

#More confidence.level assignments before finalising file.
for (i in 1:nrow(final.annotation.df)) {
  if (is.na(final.annotation.df$confidence.level[i]) && 
      (!is.na(final.annotation.df$ms2query.NPC.pathway[i]) || !is.na(final.annotation.df$canopus.NPC.pathway[i])) &&
      (!is.na(final.annotation.df$ms2query.NPC.pathway[i]) && final.annotation.df$ms2query.NPC.pathway[i] != "None")) { # Check if NOT NA before comparing
    final.annotation.df$confidence.level[i] <- "4"
  }
}

for (i in 1:nrow(final.annotation.df)) {
  if (is.na(final.annotation.df$confidence.level[i])) {
    final.annotation.df$confidence.level[i] <- "5"
  }
}

write.csv(final.annotation.df, "final-annotation-df.csv")