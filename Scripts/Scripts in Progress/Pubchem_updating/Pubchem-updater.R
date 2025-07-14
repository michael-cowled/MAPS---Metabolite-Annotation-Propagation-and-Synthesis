# --- Required Packages ---
if (!requireNamespace("rvest", quietly = TRUE)) install.packages("rvest")
if (!requireNamespace("jsonlite", quietly = TRUE)) install.packages("jsonlite")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("xml2", quietly = TRUE)) install.packages("xml2")

library(rvest)
library(jsonlite)
library(dplyr)
library(xml2)

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
    if (!is.null(response$IdentifierList$CID)) return(response$IdentifierList$CID[1]) else return(NA)
    
  } else if (type == "smiles" && property == "cids") {
    url <- paste0(base_url, "/compound/smiles/", URLencode(query), "/cids/JSON")
    response <- tryCatch(fromJSON(url), error = function(e) return(NULL))
    if (!is.null(response$IdentifierList$CID)) return(response$IdentifierList$CID[1]) else return(NA)
    
  } else if (type == "synonym" && property == "cids") {
    url <- paste0(base_url, "/compound/name/", URLencode(query), "/synonyms/JSON")
    response <- tryCatch(fromJSON(url), error = function(e) return(NULL))
    if (!is.null(response$InformationList$Information[[1]]$CID)) {
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
      list(
        SMILES = xml_text(xml_find_first(xml, ".//SMILES")),
        MolecularFormula = xml_text(xml_find_first(xml, ".//MolecularFormula")),
        MonoisotopicMass = as.numeric(xml_text(xml_find_first(xml, ".//MonoisotopicMass")))
      )
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
    Sys.sleep(0.2)
    
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

# --- Example Data ---
if (!exists("final_combined_unique_entries_filtered")) {
  final_combined_unique_entries_filtered <- data.frame(
    Best.Annotation = c("Caffeine", "Avocadyne Acetate", "Glycerol", "UnknownCompound"),
    Best.Annotation.Smiles = c(NA, NA, NA, "CCO"),
    stringsAsFactors = FALSE
  )
}

# Add columns
cols_needed <- c("CID", "Best.Annotation.Smiles", "Best.Annotation.MolecularFormula", "Best.Annotation.MonoisotopicMass")
for (col in cols_needed) {
  if (!col %in% names(final_combined_unique_entries_filtered)) {
    final_combined_unique_entries_filtered[[col]] <- if (col == "Best.Annotation.MonoisotopicMass") NA_real_ else NA_character_
  }
}

# --- Main Processing Loop ---
for (i in seq_len(nrow(final_combined_unique_entries_filtered))) {
  current_name <- final_combined_unique_entries_filtered$Best.Annotation[i]
  current_smiles <- final_combined_unique_entries_filtered$Best.Annotation.Smiles[i]
  
  if (!is.na(current_name) && current_name != "") {
    message(paste("Processing:", current_name))
    
    info <- get_cid_with_fallbacks(current_name, current_smiles)
    Sys.sleep(0.2)
    
    if (!is.null(info)) {
      final_combined_unique_entries_filtered$CID[i] <- info$CID
      if (!is.na(info$ResolvedName) && info$ResolvedName != current_name) {
        info$ResolvedName <- sub(" \\|.*", "", info$ResolvedName)
        final_combined_unique_entries_filtered$Best.Annotation[i] <- info$ResolvedName
        message(paste("  Updated Best.Annotation:", current_name, "→", info$ResolvedName))
      }
      final_combined_unique_entries_filtered$Best.Annotation.Smiles[i] <- info$SMILES
      final_combined_unique_entries_filtered$Best.Annotation.MolecularFormula[i] <- info$MolecularFormula
      final_combined_unique_entries_filtered$Best.Annotation.MonoisotopicMass[i] <- info$MonoisotopicMass
    } else {
      message(paste("  No CID found for", current_name))
    }
  }
}

#Some final tidying up
cid_cache_df$ResolvedName <- sub(" \\|.*", "", cid_cache_df$ResolvedName) #Final tidy up.
cid_cache_df <- cid_cache_df %>%
  group_by(CID) %>%
  mutate(
    ResolvedName = ifelse(is.na(ResolvedName), ResolvedName[!is.na(ResolvedName)][1], ResolvedName),
    SMILES = ifelse(is.na(SMILES), SMILES[!is.na(SMILES)][1], SMILES),
    MolecularFormula = ifelse(is.na(MolecularFormula), MolecularFormula[!is.na(MolecularFormula)][1], MolecularFormula),
    MonoisotopicMass = ifelse(is.na(MonoisotopicMass), MonoisotopicMass[!is.na(MonoisotopicMass)][1], MonoisotopicMass)
  ) %>%
  ungroup()

# Save updated cache
write.csv(cid_cache_df, cid_cache_file, row.names = FALSE)
