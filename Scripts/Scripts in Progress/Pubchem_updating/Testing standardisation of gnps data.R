###TESTING GNPS data, tidying, normalisation and probability calculations

# Applying pubchem standardisation
# --- Required Packages ---
if (!requireNamespace("rvest", quietly = TRUE)) install.packages("rvest")
if (!requireNamespace("jsonlite", quietly = TRUE)) install.packages("jsonlite")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("xml2", quietly = TRUE)) install.packages("xml2")

library(rvest)
library(jsonlite)
library(dplyr)
library(xml2)
library(dplyr)

#gnps.task.id
dataset.id <- "HGMD_0127"

excel_file <- "C:/Users/mcowled/The University of Melbourne/MA Human Gut Metabolome - Documents/Data Management and Analysis/HGM - Data Management System.xlsx"

## Updating Metadata                                                              ----> Remove from published version

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

# GNPS2 Task ID: See URL for job: https://gnps2.org/status?task=TASKID 
dataset <- read.csv("HGM/D - Dataset.csv") %>%
  filter(HGMD.ID == dataset.id)
gnps.task.id <- dataset$gnps.task.ID[1] #TASKID = unique set of random numbers and letters

#Load in gnps data
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
gnps.annotation.data$gnps.in.silico.bile.acid.info <- NA
gnps.annotation.data$gnps.in.silico.bile.acid.info[grepl("Candidate ", gnps.annotation.data$gnps.compound.name)] <- 
  gnps.annotation.data$gnps.compound.name[grepl("Candidate ", gnps.annotation.data$gnps.compound.name)]

# Apply your fix_compound_names function
if (exists("fix_compound_names") && is.function(fix_compound_names)) {
  gnps.annotation.data <- fix_compound_names(gnps.annotation.data, "gnps.compound.name")
} else {
  message("Warning: 'fix_compound_names' function not found or not a function. Skipping compound name fix.")
}

#Pre-standardisation filtering step to remove duplicates

gnps.annotation.data <- gnps.annotation.data %>%
  group_split(feature.ID) %>%
  lapply(function(df) {
    # Step 1: Filter for cosine score >= 0.7
    high_score <- df %>% filter(gnps.cosine.score >= 0.7)
    
    # If no rows have cosine score >= 0.7, keep the one with highest score
    if (nrow(high_score) == 0) {
      df <- df %>% slice_max(gnps.cosine.score, n = 1, with_ties = FALSE)
    } else {
      df <- high_score
    }
    
    # Step 2: Keep unique gnps.compound.name, then unique gnps.smiles
    df %>%
      distinct(gnps.compound.name, .keep_all = TRUE) %>%
      distinct(gnps.smiles, .keep_all = TRUE)
  }) %>%
  bind_rows()

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

# --- PubChem Query Function ---
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


# --- CID + Property Lookup with Cache Fallback ---
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
    
    new_row <- data.frame(
      LookupName = name,
      ResolvedName = ifelse(!is.na(title), title, name),
      SMILES = ifelse(!is.null(props) && !is.null(props$SMILES), props$SMILES, ifelse(!is.na(smiles), smiles, NA)),
      CID = cid,
      MolecularFormula = ifelse(!is.null(props), props$MolecularFormula, NA),
      MonoisotopicMass = ifelse(!is.null(props), props$MonoisotopicMass, NA),
      stringsAsFactors = FALSE
    )
    
    cid_cache_df <<- bind_rows(cid_cache_df, new_row)
    return(as.list(new_row[1, ]))
  }
  
  return(NULL)
}

# Add columns
cols_needed <- c("CID", "gnps.smiles", "gnps.compound.MolecularFormula", "gnps.compound.MonoisotopicMass")
for (col in cols_needed) {
  if (!col %in% names(gnps.annotation.data)) {
    gnps.annotation.data[[col]] <- if (col == "gnps.compound.MonoisotopicMass") NA_real_ else NA_character_
  }
}

# --- Main Processing Loop ---
pb <- txtProgressBar(min = 0, max = nrow(gnps.annotation.data), style = 3)

for (i in seq_len(nrow(gnps.annotation.data))) {
  current_name <- gnps.annotation.data$gnps.compound.name[i]
  current_smiles <- gnps.annotation.data$gnps.smiles[i]
  
  if (!is.na(current_name) && nzchar(current_name)) {
    message(paste("Processing:", current_name))
    
    info <- get_cid_with_fallbacks(current_name, current_smiles)
    
    if (!is.null(info)) {
      gnps.annotation.data$CID[i] <- info$CID
      if (!is.na(info$ResolvedName) && info$ResolvedName != current_name) {
        info$ResolvedName <- sub(" \\|.*", "", info$ResolvedName)
        gnps.annotation.data$gnps.compound.name[i] <- info$ResolvedName
        message(paste("  Updated gnps.compound.name:", current_name, "→", info$ResolvedName))
      }
      gnps.annotation.data$gnps.smiles[i] <- info$SMILES
      gnps.annotation.data$gnps.compound.MolecularFormula[i] <- info$MolecularFormula
      gnps.annotation.data$gnps.compound.MonoisotopicMass[i] <- info$MonoisotopicMass
    } else {
      message(paste("  No CID found for", current_name))
    }
  }
  
  # Update progress bar
  setTxtProgressBar(pb, i)
}

# Close the progress bar
close(pb)


# --- Populate from cache ---
for (i in seq_len(nrow(gnps.annotation.data))) {
  current_name <- gnps.annotation.data$gnps.compound.name[i]
  row_match <- cid_cache_df[cid_cache_df$LookupName == current_name, ]
  
  if (nrow(row_match) > 0) {
    gnps.annotation.data$CID[i] <- row_match$CID[1]
    gnps.annotation.data$gnps.smiles[i] <- row_match$SMILES[1]
    gnps.annotation.data$gnps.compound.MolecularFormula[i] <- row_match$MolecularFormula[1]
    gnps.annotation.data$gnps.compound.MonoisotopicMass[i] <- row_match$MonoisotopicMass[1]
  }
}

# --- Final Output ---
print(gnps.annotation.data)

# Save updated cache
write.csv(cid_cache_df, cid_cache_file, row.names = FALSE)

#------#
##Once standardisation has been perform, filter for unique for each feature.ID sequentially
###PER FEATURE ID SEPARATELY

#Post-standardisation filtering step to remove duplicates

gnps.annotation.data <- gnps.annotation.data %>%
  group_split(feature.ID) %>%
  lapply(function(df) {
    # Step 1: Filter for cosine score >= 0.7
    high_score <- df %>% filter(gnps.cosine.score >= 0.7)
    
    # If no rows have cosine score >= 0.7, keep the one with highest score
    if (nrow(high_score) == 0) {
      df <- df %>% slice_max(gnps.cosine.score, n = 1, with_ties = FALSE)
    } else {
      df <- high_score
    }
    
    # Step 2: Keep unique gnps.compound.name, then unique gnps.smiles
    df %>%
      distinct(gnps.compound.name, .keep_all = TRUE) %>%
      distinct(gnps.smiles, .keep_all = TRUE)
  }) %>%
  bind_rows()

##Then calculate probability for each feature with and without gold filter (gnps.library.quality)

gnps.annotation.data <- gnps.annotation.data %>%
  mutate(gnps.ID.prob = 0)

# Step 2: Assign gnps.ID.prob based on cosine score filtering
gnps.annotation.data <- gnps.annotation.data %>%
  group_by(feature.ID) %>%
  mutate(
    n_above_thresh = sum(gnps.cosine.score >= 0.7),
    gnps.ID.prob = ifelse(gnps.cosine.score >= 0.7 & n_above_thresh > 0, 1 / n_above_thresh, 0)
  ) %>%
  ungroup() %>%
  select(-n_above_thresh)  # Remove helper column

# Step 3: Filter to retain only the highest cosine score per feature.ID
gnps.annotation.data <- gnps.annotation.data %>%
  group_by(feature.ID) %>%
  arrange(desc(gnps.cosine.score)) %>%
  slice(1) %>%
  ungroup()

# OPTIONAL: Filter for distinct compounds, and return the number
gnps.annotation.data <- gnps.annotation.data %>%
  distinct(gnps.compound.name, .keep_all = TRUE)
nrow(gnps.annotation.data)