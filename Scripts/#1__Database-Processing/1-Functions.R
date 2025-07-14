##Annotation-table: Functions and Required Packages

#Removed from main script for simplicity

#------------------------------------------------------------------------------#
## 1. check_and_install
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

#------------------------------------------------------------------------------#
## 2. fix_compound_names 
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

#------------------------------------------------------------------------------#
## 3. get_pubchem
get_pubchem <- function(query, type, property = NULL) {
  base_url <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
  
  if (type == "name" && property == "cids") {
    url <- paste0(base_url, "/compound/name/", URLencode(query), "/cids/JSON")
    response <- tryCatch(fromJSON(url), error = function(e) return(NULL))
    if (!is.null(response$IdentifierList$CID)) {
      Sys.sleep(0.2)
      return(response$IdentifierList$CID[1])
    } else return(NA)
    
  } else if (type == "smiles" && property == "cids") {
    url <- paste0(base_url, "/compound/smiles/", URLencode(query), "/cids/JSON")
    response <- tryCatch(fromJSON(url), error = function(e) return(NULL))
    if (!is.null(response$IdentifierList$CID)) {
      Sys.sleep(0.2)
      return(response$IdentifierList$CID[1])
    } else return(NA)
    
  } else if (type == "synonym" && property == "cids") {
    url <- paste0(base_url, "/compound/name/", URLencode(query), "/synonyms/JSON")
    response <- tryCatch(fromJSON(url), error = function(e) return(NULL))
    if (!is.null(response$InformationList$Information[[1]]$CID)) {
      Sys.sleep(0.2)
      return(response$InformationList$Information[[1]]$CID)
    } else {
      message("    Synonym search returned no CID.")
      return(NA)
    }
    
  } else if (type == "cid" && property == "title") {
    page_url <- paste0("https://pubchem.ncbi.nlm.nih.gov/compound/", query)
    page <- tryCatch(read_html(page_url), error = function(e) return(NULL))
    if (!is.null(page)) {
      title <- html_text(html_node(page, "title"))
      return(gsub(" - PubChem", "", title, fixed = TRUE))
    } else return(NA)
    
  } else if (type == "cid" && property == "properties") {
    url <- paste0(base_url, "/compound/cid/", query, "/property/SMILES,MolecularFormula,MonoisotopicMass/XML")
    response <- tryCatch({
      xml <- read_xml(url)
      result <- list(
        SMILES = xml_text(xml_find_first(xml, ".//SMILES")),
        MolecularFormula = xml_text(xml_find_first(xml, ".//MolecularFormula")),
        MonoisotopicMass = as.numeric(xml_text(xml_find_first(xml, ".//MonoisotopicMass")))
      )
      Sys.sleep(0.2)
      return(result)
    }, error = function(e) return(NULL))
    return(response)
    
  } else {
    message("Invalid query type or property.")
    return(NA)
  }
}

#------------------------------------------------------------------------------#
## 4. get_cid_with_fallbacks
get_cid_with_fallbacks <- function(name, smiles = NA) {
  # Prefer to match by exact name first
  name_match <- cid_cache_df[!is.na(cid_cache_df$LookupName) & cid_cache_df$LookupName == name, ]
  
  # If no match on name, try matching on SMILES (if available)
  if (nrow(name_match) == 0 && !is.na(smiles) && smiles != "") {
    name_match <- cid_cache_df[!is.na(cid_cache_df$SMILES) & cid_cache_df$SMILES == smiles, ]
  }
  
  if (nrow(name_match) > 0 && !is.na(name_match$CID[1])) {
    message(paste("  [CACHE] CID found for:", name))
    return(list(
      CID = name_match$CID[1],
      ResolvedName = name_match$ResolvedName[1],
      SMILES = name_match$SMILES[1],
      MolecularFormula = name_match$MolecularFormula[1],
      MonoisotopicMass = name_match$MonoisotopicMass[1]
    ))
  }
  
  # Step 1: Name → CID
  cid <- get_pubchem(name, "name", "cids")
  
  # Step 2: Try SMILES
  if (is.na(cid) && !is.na(smiles) && smiles != "") {
    message(paste("  Name lookup failed. Trying SMILES:", smiles))
    cid <- get_pubchem(smiles, "smiles", "cids")
  }
  
  # Step 3: Try synonym
  if (is.na(cid)) {
    message(paste("  Name and SMILES failed. Trying synonym search for:", name))
    cid <- get_pubchem(name, "synonym", "cids")
  }
  
  # If CID found, get details and cache
  if (!is.na(cid)) {
    props <- get_pubchem(cid, "cid", "properties")
    title <- get_pubchem(cid, "cid", "title")
    
    # Ensure props is a list before accessing
    if (!is.list(props)) props <- list()
    if (!is.character(title)) title <- NA
    
    new_row <- data.frame(
      LookupName = name,
      ResolvedName = if (!is.null(title) && !is.na(title)) title else name,
      SMILES = if (!is.null(props$SMILES)) props$SMILES else if (!is.na(smiles)) smiles else NA,
      CID = cid,
      MolecularFormula = if (!is.null(props$MolecularFormula)) props$MolecularFormula else NA,
      MonoisotopicMass = if (!is.null(props$MonoisotopicMass)) props$MonoisotopicMass else NA,
      stringsAsFactors = FALSE
    )
    
    cid_cache_df <<- bind_rows(cid_cache_df, new_row)
    return(as.list(new_row[1, ]))
  }
  
  return(NULL)
}

#------------------------------------------------------------------------------#
## 5. standardise_annotation
standardise_annotation <- function(data, name_col, smiles_col) {
  # Filter out rows where name_col contains "candidate" (case-insensitive)
  data <- data[!grepl("candidate", data[[name_col]], ignore.case = TRUE), ]
  
  n <- nrow(data)
  pb <- txtProgressBar(min = 0, max = n, style = 3)
  
  # Ensure required output columns exist
  data[[paste0("CID")]] <- NA
  data[[paste0("MolecularFormula")]] <- NA
  data[[paste0("MonoisotopicMass")]] <- NA
  
  for (i in seq_len(n)) {
    current_name <- data[[name_col]][i]
    current_smiles <- data[[smiles_col]][i]
    
    if (!is.na(current_name) && nzchar(current_name)) {
      message(paste("Processing:", current_name))
      
      info <- get_cid_with_fallbacks(current_name, current_smiles)
      
      if (!is.null(info)) {
        data[[paste0("CID")]][i] <- info$CID
        if (!is.na(info$ResolvedName) && info$ResolvedName != current_name) {
          info$ResolvedName <- sub(" \\|.*", "", info$ResolvedName)
          data[[name_col]][i] <- info$ResolvedName
          message(paste("  Updated", name_col, ":", current_name, "→", info$ResolvedName))
        }
        data[[smiles_col]][i] <- info$SMILES
        data[[paste0("MolecularFormula")]][i] <- info$MolecularFormula
        data[[paste0("MonoisotopicMass")]][i] <- info$MonoisotopicMass
      } else {
        message(paste("  No CID found for", current_name))
      }
    }
    
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  return(data)
}

#------------------------------------------------------------------------------#
## 6. get_result
# Function to get the first non-NA value from specified columns, the column name, and the superclass
get_result <- function(paired_value, summary_data) {
  # Define the compound name column
  compound_col <- "compound.name"
  
  # Define the superclass columns to check
  superclass_columns <- c("NPC.superclass", "canopus.NPC.superclass")
  
  # Filter and ensure only a single row is returned
  data_subset <- summary_data %>%
    filter(feature.ID == paired_value) %>%
    select(all_of(compound_col), all_of(superclass_columns), smiles) %>%
    slice(1)  # Take the first row if duplicates exist
  
  compound_value <- data_subset[[compound_col]]
  
  if (!is.na(compound_value) && compound_value != "null") {
    # Find the first valid superclass
    superclass_value <- NA
    for (col in superclass_columns) {
      sc_val <- data_subset[[col]]
      if (!is.na(sc_val) && sc_val != "null") {
        superclass_value <- sc_val
        break
      }
    }
    return(list(value = compound_value, column = compound_col, superclass = superclass_value))
  }
  
  # Fallback to smiles
  csi_value <- data_subset$smiles
  if (!is.na(csi_value) && csi_value != "null") {
    return(list(value = csi_value, column = "smiles", superclass = NA))
  }
  
  return(list(value = NA, column = NA, superclass = NA))
}

#------------------------------------------------------------------------------#  
## 7. paired_feature_finder
# Function: Finds the features linked by an edge in a GNPS cluster
paired_feature_finder <- function(ID) {
  filtered_pairs <- gnps.cluster.pairs %>%
    filter(CLUSTERID1 == ID | CLUSTERID2 == ID) %>%
    arrange(desc(Cosine)) # Sort by decreasing Cosine
  
  # Extract paired values
  paired_values <- filtered_pairs %>%
    mutate(paired_value = ifelse(CLUSTERID1 == ID, CLUSTERID2, CLUSTERID1)) %>%
    pull(paired_value)
  
  return(paired_values)
}

#------------------------------------------------------------------------------#  
# 8. process_iin_group
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

#------------------------------------------------------------------------------#  
## 9. update_data_frame
# Define a function to update the main data frame
update_data_frame <- function(df, results, id_column) {
  df_non_na <- df %>% filter(!is.na(!!sym(id_column)))
  results_df <- bind_rows(results)
  
  updated_df <- df_non_na %>% filter(feature.ID %in% results_df$feature.ID)
  df_na <- df %>% filter(is.na(!!sym(id_column)))
  
  final_df <- bind_rows(updated_df, df_na)
  return(final_df)
}


## 10. process_all_features
# Main processing function for collapsing by ion identity networking
process_all_features <- function(data, id_column, feature_columns) {
  unique_ids <- unique(data[[id_column]][!is.na(data[[id_column]])])
  results <- lapply(unique_ids, function(id) {
    group_data <- filter(data, !!sym(id_column) == id)
    process_iin_group(group_data, feature_columns)
  })
  names(results) <- as.character(unique_ids)
  return(results)
}
#------------------------------------------------------------------------------#  
## 11. process_data_fixed
redundancy_fixer <- function(dataset, column_to_check = "Best.Annotation.Smiles", rt_column = "rt", rt_tolerance = 1) {
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

#------------------------------------------------------------------------------#  
## 12. validate_and_get_paths
validate_and_get_paths <- function(folder) {
  # Clean up the folder path (replace backslashes with forward slashes)
  # Note the double backslashes in the pattern for escaping the backslash character itself
  folder <- gsub("\\\\", "/", folder)
  
  if (dir.exists(folder)) {
  } else {
    stop(paste0("The folder", folder, "does not exist."))  # Stop execution
  }
  
  # Create a named list of paths for easier iteration and return
  file_paths <- list(
    mzmine_data = mzmine.data,
    mzmine_annotations = mzmine.annotations,
    canopus_data = canopus.data,
    csi_data = csi.data,
    zodiac_data = zodiac.data,
    ms2query_data = ms2query.data,
    cytoscape = cytoscape
  )
  
  # Define specific error messages for each file
  error_messages <- list(
    mzmine_data = "The file or folder is missing for mzmine data (ms1-and-ms2.csv).",
    mzmine_annotations = "The file or folder is missing for mzmine data (data_annotations.csv).",
    canopus_data = "The file or folder is missing for canopus data (canopus_structure_summary.tsv).",
    csi_data = "The file or folder is missing for csi:fingerID data (structure_identifications_top-100.tsv). Recompute for top K=100 hits.",
    zodiac_data = "The file or folder is missing for zodiac data (formula_identifications.tsv).",
    ms2query_data = "The file or folder is missing for ms2query data (ms2query.csv).",
    cytoscape = "The file or folder is missing for cytoscape.csv."
  )
  
  # Check if each file exists
  for (name in names(file_paths)) {
    if (!file.exists(file_paths[[name]])) {
      stop(error_messages[[name]]) # Stop execution with the specific error message
    }
  }
  
  # If all checks pass, return the list of file paths
  return(file_paths)
}

#------------------------------------------------------------------------------#  
## 13. filter_data_by_standards
filter_data_by_standards <- function(
    data,              # either gnps.data or csi.data
    mzmine.data,       # standards data frame (same as before)
    compound_col,      # column name for compound ID or CID in 'data'
    feature_id_col = "feature.ID",
    rt_tolerance = 0.1
) {
  library(dplyr)
  
  # Step 1: Filter authentic standards from mzmine.data
  mzmine.standards <- filter(mzmine.data, !is.na(compound.name))
  
  # Step 2: Unique CIDs to check (from standards)
  mzmine.to.check <- unique(mzmine.standards$CID)
  
  # Step 3: Join RT from mzmine.data, ensure numeric
  data.cleaned <- data %>%
    left_join(mzmine.data %>% select(all_of(feature_id_col), rt), by = feature_id_col) %>%
    mutate(rt = suppressWarnings(as.numeric(as.character(rt))))
  
  # Step 4: Loop through each CID in standards to filter data.cleaned
  for (i in mzmine.to.check) {
    temp <- filter(data.cleaned, .data[[compound_col]] == i)
    if (nrow(temp) == 0) next
    
    temp2 <- filter(mzmine.standards, CID == i)
    temp3 <- filter(temp, .data[[feature_id_col]] != temp2$feature.ID[1])  # standard feature.ID
    
    rt.std.raw <- temp2$rt[1]
    rt.std <- suppressWarnings(as.numeric(as.character(rt.std.raw)))
    
    if (length(rt.std) == 0 || is.na(rt.std)) {
      warning(paste("Skipping CID", i, "- invalid standard rt:", rt.std.raw))
      next
    }
    
    # Remove entries with RT difference > rt_tolerance
    for (j in temp3[[feature_id_col]]) {
      rt.j.raw <- data.cleaned$rt[data.cleaned[[feature_id_col]] == j]
      rt.j <- suppressWarnings(as.numeric(as.character(rt.j.raw)))
      if (length(rt.j) == 0 || all(is.na(rt.j))) next
      
      if (abs(rt.std - rt.j[1]) > rt_tolerance) {
        data.cleaned <- filter(data.cleaned, .data[[feature_id_col]] != j)
      }
    }
    
    # Remove remaining entries within tolerance but not standard
    temp <- filter(data.cleaned, .data[[compound_col]] == i)
    temp3 <- filter(temp, .data[[feature_id_col]] != temp2$feature.ID[1])
    
    rt.others <- vapply(temp3[[feature_id_col]], function(fid) {
      rt.val <- data.cleaned$rt[data.cleaned[[feature_id_col]] == fid]
      if (length(rt.val) == 0 || all(is.na(rt.val))) return(NA_real_)
      suppressWarnings(as.numeric(as.character(rt.val[1])))
    }, numeric(1))
    
    valid_idx <- !is.na(rt.others)
    
    close_match_ids <- temp3[[feature_id_col]][valid_idx & abs(rt.std - rt.others[valid_idx]) <= rt_tolerance]
    
    data.cleaned <- filter(data.cleaned, !(.data[[feature_id_col]] %in% close_match_ids))
  }
  
  return(data.cleaned)
}

#------------------------------------------------------------------------------#  
## 14. deduplicate_data
# Default values for gnps.data (but should work for others)
deduplicate_data <- function(data, name_col, score_col) {
  data %>%
    arrange(desc({{score_col}})) %>%
    group_by({{name_col}}) %>%
    slice_head(n = 1) %>%
    ungroup()
}

