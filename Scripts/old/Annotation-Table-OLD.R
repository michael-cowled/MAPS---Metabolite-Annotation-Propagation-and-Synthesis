#### ANNOTATION TABLE SCRIPT FOR AUSTRALIAN HUMAN GUT METABOLOME DATABASE ####

#------------------------------------------------------------------------------#
              ###---START OF USER-FED INFORMATION---###

# Dataset ID: From Data Management Plan:

  dataset.id <- "HGMD_0110"     #####Change to dataset of interest#####

# Specify the path to the Data Management Plan (unless moved to Mediaflux then should be consistent)

  excel_file <- "C:/Users/mcowled/The University of Melbourne/MA Human Gut Metabolome - Documents/Data Management and Analysis/HGM - Data Management System.xlsx"

# Annotation Acceptance Probabilities (Defaults for Exploratory Analyses)
# Increase stringency if the application demands it.

  canopus.prob <- 0.7           # Default is 0.7
  csi.prob <- 0.8863            # Default is 0.8863
  ms2query.prob <- 0.63         # Default is 0.63

              ###---END OF USER-FED INFORMATION---###
#------------------------------------------------------------------------------#
  
  
## Loading libraries

  # Function to check and install missing packages, then load them
  check_and_install <- function(packages) {
    to_install <- packages[!(packages %in% installed.packages()[, "Package"])]
    if (length(to_install) > 0) {
      message(paste("Installing:", paste(to_install, collapse = ", ")))
      install.packages(to_install, dependencies = TRUE)
      library(packages)
    }
    # Always load the packages, regardless of whether they were just installed
    message(paste("Loading:", paste(packages, collapse = ", ")))
    invisible(lapply(packages, library, character.only = TRUE))
  }
  
  required_packages <- c("dplyr", "tidyr", "stringr", "readr", "reshape2", "ggplot2", "svglite", "readxl", "openxlsx", "tidyverse")
  check_and_install(required_packages)

  
  
## Updating Metadata                                                              ----> Remove from published version

  sheet_names <- excel_sheets(excel_file) # Get the sheet names
  for (sheet in sheet_names) {
    data <- read_excel(excel_file, sheet = sheet)
    write.csv(data, file = paste0("HGM/", sheet, ".csv"), row.names = FALSE)
  }



# GNPS2 Task ID: See URL for job: https://gnps2.org/status?task=TASKID 
dataset <- read.csv("HGM/D - Dataset.csv") %>%
  filter(HGMD.ID == dataset.id)
gnps.task.id <- dataset$gnps.task.ID[1] #TASKID = unique set of random numbers and letters

if (is.na(gnps.task.id)) {
  stop("gnps.task.ID is missing. Add to 'dataset' csv before re-running.")  # Stop execution
}

#File List Extractor#
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

## Check if all data is present in folder with correct naming conventions         ###CHECK SPELLING!!!!###
#Folder with mzmine, ms2query, sirius and gnps results
folder <- dataset$Processed.Data.Folder[1]  ##Determines folder to process based on dataset.id

if (dir.exists(folder)) {
} else {
  stop(paste0("The folder", folder, "does not exist."))  # Stop execution
}
# Clean up the folder path (replace backslashes with forward slashes):
folder <- gsub("\\\\", "/", folder) # Note the double backslashes in the pattern

## Import annotation tables to merge # mzmine - works for v4.5.0

  mzmine.data <- paste0(folder, "/mzmine/ms1-and-ms2.csv")
  canopus.data <- paste0(folder, "/sirius/canopus_structure_summary.tsv")
  csi.data <- paste0(folder, "/sirius/structure_identifications.tsv")
  zodiac.data <- paste0(folder, "/sirius/formula_identifications.tsv")
  ms2query.data <- paste0(folder, "/ms2query/ms2query.csv")
  cytoscape <- paste0(folder, "/gnps/cytoscape.csv")
  
  if (file.exists(mzmine.data)) {
  } else {
    stop(paste0("The file or folder is missing for mzmine data."))  # Stop execution
  }
  
  if (file.exists(canopus.data)) {
  } else {
    stop(paste0("The file or folder is missing for canopus data."))  # Stop execution
  }
  
  if (file.exists(csi.data)) {
  } else {
    stop(paste0("The file or folder is missing for csi:fingerID data."))  # Stop execution
  }
  
  if (file.exists(zodiac.data)) {
  } else {
    stop(paste0("The file or folder is missing for zodiac data."))  # Stop execution
  }
  
  if (file.exists(ms2query.data)) {
  } else {
    stop(paste0("The file or folder is missing for ms2query data."))  # Stop execution
  }
  
  if (file.exists(cytoscape)) {
  } else {
    stop(paste0("The file or folder is missing for cytoscape.csv."))  # Stop execution
  }

##Load in MZMINE data
mzmine.data <- read.csv(mzmine.data) # Derived from Export to CSV file (modular)
sample.data <- mzmine.data # A copy to be used for different processing later
if (!("spectral_db_matches.compound_name" %in% names(mzmine.data))) {
  mzmine.data$spectral_db_matches.compound_name <- NA
}
mzmine.data <- select(mzmine.data, "id", "rt", "mz", "ion_identities.iin_id", "spectral_db_matches.compound_name") %>%
  mutate(across(everything(), as.character)) %>%
  mutate(across(everything(), ~ na_if(., "")))
names(mzmine.data) <- c('feature.ID', "rt", "mz", "ion.identity.ID", "authentic.standard")
mzmine.data$feature.ID <- as.numeric(mzmine.data$feature.ID)
mzmine.data$authentic.standard.smiles <- NA                         

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
        
        fix_compound_names <- function(df, column) {
          # Convert to data.table for faster operations
          df <- data.table::as.data.table(df)
          
          # Define patterns and replacements for initial cleaning
          patterns <- c(
            "\\.(alpha|beta|gamma|omega)\\.-", "\\.(DELTA)\\.", "[[:space:]]*\\[IIN-based(?:[[:space:]]+on)?[[:space:]]*:[[:space:]]*[^)]*\\]",
            "NCG[[:alnum:]-]+_[[:alnum:]]+_", "NCG[^\\s!]+[\\s!]?",
            "Spectral Match to ", " from NIST14",
            "40.0 eV", "70.0 eV", "50.0 eV", ".00 eV", '"', "CollisionEnergy", "\\|.*",
            "_", "\\[", "\\]", "\\{", "\\}", " :102040", " :205060", " - ",
            "30.0 eV", "Massbank:[[:alnum:]]+\\s? ", "MassbankEU:[[:alnum:]]+\\s? ",
            "MoNA:[[:alnum:]]+\\s? ", " (Chimeric precursor selection)", "; \\(M\\+.*"
          )
          replacements <- c(
            "\\1-", "delta-", "", "", "", "", "", "", "", "", "", "", "", "",
            " ", "(", ")", "(", ")", "", "", "", "", "", "", "", "", "", ""
          )
          
          # Loop through columns and apply initial fixes
          for (col in column) {
            for (i in seq_along(patterns)) {
              df[grepl(patterns[i], get(col), ignore.case = TRUE), (col) := gsub(patterns[i], replacements[i], get(col))]
            }
            
            # --- New section for capitalization rules ---
            # Apply capitalization rules after initial cleaning
            
            # 1. If ALL CAPS, change to first letter capital, rest lowercase
            # Check if the entire string is uppercase (excluding spaces and non-alphabetic characters for the check)
            df[, (col) := ifelse(
              grepl("^[[:upper:]\\s\\d\\W]*$", get(col)) & grepl("[[:alpha:]]", get(col)), # Ensure it contains at least one letter
              paste0(toupper(substring(get(col), 1, 1)), tolower(substring(get(col), 2))),
              get(col)
            )]
            
            # 2. If ALL LOWERCASE, change first letter to uppercase
            # Check if the entire string is lowercase (excluding spaces and non-alphabetic characters for the check)
            df[, (col) := ifelse(
              grepl("^[[:lower:]\\s\\d\\W]*$", get(col)) & grepl("[[:alpha:]]", get(col)), # Ensure it contains at least one letter
              paste0(toupper(substring(get(col), 1, 1)), substring(get(col), 2)),
              get(col)
            )]
            
            # Mixed case names are left unchanged by these rules, as they won't match the ALL CAPS or ALL LOWERCASE patterns.
            # --- End of new section ---
          }
          
          return(as.data.frame(df))
        }
        
        # Apply the function to your data frame
        gnps.data <- fix_compound_names(gnps.data, "gnps.compound.name")
        
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
csi.data$csi.compound.name[grepl("Solaparnaine", csi.data$csi.compound.name, ignore.case = TRUE)] <- "Solaparnaine" ##troublesome case

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

# MS2QUERY:write
ms2query.data <- read.csv(ms2query.data)
ms2query.data  <- ms2query.data[, c(2, 3, 7, 8, 10, 17, 18)]
names(ms2query.data) <- c("ms2query.score", "ms2query.mzdiff", 
                          "ms2query.analogue.compound.name", "ms2query.smiles", 
                          'feature.ID', "ms2query.NPC.superclass", "ms2query.NPC.pathway")
ms2query.data <- ms2query.data[, c(3, 1, 2, 4, 5, 7, 6)] #Reorders columns
ms2query.data <- fix_compound_names(ms2query.data, "ms2query.analogue.compound.name")  ##df, and column to fix


## Merge annotations into one big table                                     
full.annotation.data <- mzmine.data %>%
  full_join(gnps.data, by = "feature.ID") %>%
  full_join(canopus.data, by = "feature.ID") %>%
  full_join(csi.data, by = "feature.ID") %>%
  full_join(zodiac.data, by = "feature.ID") %>%
  full_join(ms2query.data, by = "feature.ID")

#PubChem Scraper - Finds SMILES if missing

library(httr)

# --- Configuration ---
SMILES_CACHE_FILE <- "smiles_cache.csv" # Name of your cache file

# --- Function to load or create SMILES cache ---
load_smiles_cache <- function(file_path) {
  if (file.exists(file_path)) {
    message("Loading existing SMILES cache from: ", file_path)
    return(read.csv(file_path, stringsAsFactors = FALSE))
  } else {
    message("Creating new SMILES cache at: ", file_path)
    return(data.frame(
      compound_name = character(),
      smiles = character(),
      stringsAsFactors = FALSE
    ))
  }
}

# --- Function to save SMILES cache ---
save_smiles_cache <- function(cache_df, file_path) {
  write.csv(cache_df, file_path, row.names = FALSE)
  message("SMILES cache saved to: ", file_path)
}

# --- Function to fetch Synonyms from PubChem, handling errors ---
get_pubchem <- function(name, input, output) { # Typical inputs = 'name', 'cid'; outputs = 'property/SMILES', 'synonyms', 'cids'
  if (is.na(name) || name == "") { # Add check for empty string
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
      content <- content(response, "text", encoding = "UTF-8") # Read content as text
      # PubChem API sometimes returns "No CID Found" or similar messages
      if (grepl("No CID Found|Not Found", content, ignore.case = TRUE)) {
        return(NA)
      }
      return(as.character(gsub("\n.*", "", content))) # Clean up multi-line responses
    }
  }, error = function(e) {
    warning(paste("Error fetching data for", name, ":", e$message))
    return(NA)
  })
}

# --- Main Script Execution ---

# 1. Load or create the SMILES cache
smiles_cache <- load_smiles_cache(SMILES_CACHE_FILE)

# --- Process 'authentic.standard' column ---
message("\n--- Processing 'authentic.standard' SMILES ---")
for (i in seq_len(nrow(full.annotation.data))) {
  compound_name <- full.annotation.data$authentic.standard[i]
  
  if (!is.na(compound_name) &&
      (is.na(full.annotation.data$authentic.standard.smiles[i]) || full.annotation.data$authentic.standard.smiles[i] == "N/A")) {
    
    # Check cache first
    cached_smiles <- smiles_cache$smiles[match(compound_name, smiles_cache$compound_name)]
    
    if (!is.na(cached_smiles) && cached_smiles != "N/A") {
      full.annotation.data$authentic.standard.smiles[i] <- cached_smiles
      message(paste0("Cached SMILES found for 'authentic.standard': ", compound_name, " -> ", cached_smiles))
    } else {
      # Scrape if not in cache or if cached as "N/A"
      message(paste0("Scraping 'authentic.standard': ", compound_name))
      cid <- get_pubchem(compound_name, "name", "cids")
      Sys.sleep(0.2) # Be polite to the API
      
      if (!is.na(cid)) {
        smiles <- get_pubchem(cid, "cid", "property/SMILES")
        Sys.sleep(0.2) # Be polite to the API
        full.annotation.data$authentic.standard.smiles[i] <- smiles
        
        # Add to cache (or update if previously "N/A")
        if (compound_name %in% smiles_cache$compound_name) {
          smiles_cache$smiles[smiles_cache$compound_name == compound_name] <- smiles
        } else {
          smiles_cache <- bind_rows(smiles_cache, data.frame(compound_name = compound_name, smiles = smiles, stringsAsFactors = FALSE))
        }
        save_smiles_cache(smiles_cache, SMILES_CACHE_FILE) # Save after each successful scrape
        message(paste0("Scraped SMILES for 'authentic.standard': ", compound_name, " -> ", smiles))
      } else {
        full.annotation.data$authentic.standard.smiles[i] <- "N/A"
        message(paste0("No CID/SMILES found for 'authentic.standard': ", compound_name))
        # Also add to cache as N/A to avoid repeated scraping attempts
        if (compound_name %in% smiles_cache$compound_name) {
          smiles_cache$smiles[smiles_cache$compound_name == compound_name] <- "N/A"
        } else {
          smiles_cache <- bind_rows(smiles_cache, data.frame(compound_name = compound_name, smiles = "N/A", stringsAsFactors = FALSE))
        }
        save_smiles_cache(smiles_cache, SMILES_CACHE_FILE)
      }
    }
  } else {
    message(paste0("Skipping 'authentic.standard' row ", i, ": already has SMILES or no name."))
  }
}

# --- Process 'gnps.compound.name' column ---
message("\n--- Processing 'gnps.compound.name' SMILES ---")
for (i in seq_len(nrow(full.annotation.data))) {
  compound_name <- full.annotation.data$gnps.compound.name[i]
  
  if (!is.na(compound_name) &&
      (is.na(full.annotation.data$gnps.smiles[i]) || full.annotation.data$gnps.smiles[i] == "N/A") &&
      (!grepl("Candidate", compound_name, ignore.case = TRUE))) {
    
    # Check cache first
    cached_smiles <- smiles_cache$smiles[match(compound_name, smiles_cache$compound_name)]
    
    if (!is.na(cached_smiles) && cached_smiles != "N/A") {
      full.annotation.data$gnps.smiles[i] <- cached_smiles
      message(paste0("Cached SMILES found for 'gnps.compound.name': ", compound_name, " -> ", cached_smiles))
    } else {
      # Scrape if not in cache or if cached as "N/A"
      message(paste0("Scraping 'gnps.compound.name': ", compound_name))
      cid <- get_pubchem(compound_name, "name", "cids")
      Sys.sleep(0.2) # Be polite to the API
      
      if (!is.na(cid)) {
        smiles <- get_pubchem(cid, "cid", "property/SMILES")
        Sys.sleep(0.2) # Be polite to the API
        full.annotation.data$gnps.smiles[i] <- smiles
        
        # Add to cache (or update if previously "N/A")
        if (compound_name %in% smiles_cache$compound_name) {
          smiles_cache$smiles[smiles_cache$compound_name == compound_name] <- smiles
        } else {
          smiles_cache <- bind_rows(smiles_cache, data.frame(compound_name = compound_name, smiles = smiles, stringsAsFactors = FALSE))
        }
        save_smiles_cache(smiles_cache, SMILES_CACHE_FILE) # Save after each successful scrape
        message(paste0("Scraped SMILES for 'gnps.compound.name': ", compound_name, " -> ", smiles))
      } else {
        full.annotation.data$gnps.smiles[i] <- "N/A"
        message(paste0("No CID/SMILES found for 'gnps.compound.name': ", compound_name))
        # Also add to cache as N/A to avoid repeated scraping attempts
        if (compound_name %in% smiles_cache$compound_name) {
          smiles_cache$smiles[smiles_cache$compound_name == compound_name] <- "N/A"
        } else {
          smiles_cache <- bind_rows(smiles_cache, data.frame(compound_name = compound_name, smiles = "N/A", stringsAsFactors = FALSE))
        }
        save_smiles_cache(smiles_cache, SMILES_CACHE_FILE)
      }
    }
  } else {
    message(paste0("Skipping 'gnps.compound.name' row ", i, ": already has SMILES, no name, or is 'Candidate'."))
  }
}

message("\n--- Scraping process complete ---")

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
    Propagated.Annotation.Type = NA, # Initialize the new column
    Propagated.Annotation.Class = NA # This is already here
  )

# Function to get the first non-NA value from specified columns, the column name, and the superclass
get_result <- function(paired_value, summary_data) {
  # Define the columns to check for non-NA values
  columns_to_check <- c("gnps.compound.name", "csi.compound.name",
                        "ms2query.analogue.compound.name")
  
  # Define the corresponding superclass columns
  superclass_columns <- c("gnps.NPC.superclass", "canopus.NPC.superclass",
                          "ms2query.NPC.superclass")
  
  # Filter the relevant row and select the columns
  data_subset <- summary_data %>%
    filter(feature.ID == paired_value) %>%
    select(all_of(columns_to_check), all_of(superclass_columns), csi.smiles)
  
  # Find the first non-NA and non-'null' value and its corresponding superclass
  result_value <- NA
  column_name <- NA
  superclass_value <- NA
  
  for (i in seq_along(columns_to_check)) {
    col <- columns_to_check[i]
    superclass_col <- superclass_columns[i]
    
    value <- data_subset[[col]]
    if (!is.na(value) && value != 'null') {
      result_value <- value
      column_name <- col
      superclass_value <- data_subset[[superclass_col]] # Get the corresponding superclass
      break
    }
  }
  
  # If no valid result is found from compound names, fallback to 'csi.smiles'
  # In this case, there's no direct superclass to propagate from csi.smiles itself.
  # You might want to decide what 'superclass_value' should be in this case (e.g., NA, "Unknown", etc.)
  if (is.na(result_value) || result_value == 'null') {
    csi_value <- data_subset$csi.smiles
    if (!is.na(csi_value) && csi_value != 'null') {
      result_value <- csi_value
      column_name <- "csi.smiles"
      superclass_value <- NA # Or "Unknown" or similar, if appropriate
    }
  }
  
  # Return the result, corresponding column name, and the superclass value
  list(value = result_value, column = column_name, superclass = superclass_value)
}

# The rest of the code remains the same:
# Function: Finds the features linked by an edge in a GNPS cluster
paired.feature.finder <- function(ID) {
  filtered_pairs <- gnps.cluster.pairs %>%
    filter(CLUSTERID1 == ID | CLUSTERID2 == ID) %>%
    arrange(desc(Cosine)) # Sort by decreasing Cosine
  
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

results <- process_all_features(data, "ion.identity.ID", c("feature1", "feature2"))

# Print results for each ID
for (id in names(results)) {
  cat("Results for ion.identity.ID =", id, ":\n")
  print(results[[id]])
}

# Update propagation.df by removing rows based on results and retaining NA IDs
propagation.df <- update_data_frame(propagation.df, results, "ion.identity.ID")

####Cytoscape Tidier####

##Uses output from Annotation-Table.R##

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

##Redundancy reduction

          dataset <- final.annotation.df %>%
            filter(!is.na(Best.Annotation.Smiles) & Best.Annotation.Smiles != "N/A") # Changed OR to AND
          dataset$rt <- as.numeric(dataset$rt)
          dataset$redundant <- FALSE
          
          other.data <- final.annotation.df %>%
            filter(is.na(Best.Annotation.Smiles) | Best.Annotation.Smiles == "N/A")
          
          process_data_fixed <- function(dataset, column_to_check = "Best.Annotation.Smiles", rt_column = "rt", rt_tolerance = 1) {
            # --- Error Handling: Check column existence ---
            if (!column_to_check %in% colnames(dataset)) {
              stop(paste("Column '", column_to_check, "' not found in the dataset.", sep = ""))
            }
            if (!rt_column %in% colnames(dataset)) {
              stop(paste("Column '", rt_column, "' not found in the dataset.", sep = ""))
            }
            
            # --- Ensure rt_column is numeric ---
            if (!is.numeric(dataset[[rt_column]])) {
              warning(paste("Column '", rt_column, "' is not numeric. Attempting to convert.", sep = ""))
              dataset[[rt_column]] <- as.numeric(as.character(dataset[[rt_column]]))
              if (any(is.na(dataset[[rt_column]]))) {
                stop(paste("Some values in ", rt_column, " could not be converted to numeric. Please inspect data", sep = ""))
              }
            }
            
            # Add a new column to mark redundancy (initialized as FALSE)
            dataset$redundant <- FALSE
            
            unique_values <- unique(dataset[[column_to_check]]) # corrected line
            
            for (value in unique_values) {
              temp_df <- dataset %>% filter(!!sym(column_to_check) == value)
              
              if (nrow(temp_df) > 1) {
                rt_values <- sort(temp_df[[rt_column]])
                rt_values <- rt_values[!is.na(rt_values)]
                
                grouped_rts <- list()
                used_indices <- rep(FALSE, length(rt_values))
                
                for (i in seq_along(rt_values)) {
                  if (!used_indices[i]) {
                    group <- rt_values[i]
                    used_indices[i] <- TRUE
                    for (j in (i + 1):length(rt_values)) {
                      if (!is.na(rt_values[i]) && !is.na(rt_values[j])) {
                        if (abs(rt_values[j] - rt_values[i]) <= rt_tolerance) {
                          group <- c(group, rt_values[j])
                          used_indices[j] <- TRUE
                        }
                      }
                    }
                    grouped_rts <- c(grouped_rts, list(sort(group)))
                  }
                }
                
                for (group in grouped_rts) {
                  group_df <- temp_df %>% filter(!!sym(rt_column) %in% group)
                  if (nrow(group_df) > 1) {
                    dataset$redundant[which(dataset[[rt_column]] %in% group_df[[rt_column]][-1])] <- TRUE
                  }
                }
              }
            }
            
            return(dataset)
          }
                  
        dataset <- process_data_fixed(dataset, column_to_check = "Best.Annotation.Smiles")
        
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