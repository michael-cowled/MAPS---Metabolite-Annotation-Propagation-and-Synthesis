# --- Load required packages ---
library(dplyr)
library(jsonlite)
library(rvest)
library(xml2)

# --- Load existing cache ---
cid_cache_file <- "cid_cache.csv"
if (file.exists(cid_cache_file)) {
  cid_cache_df <- read.csv(cid_cache_file, stringsAsFactors = FALSE)
} else {
  stop("Cache file not found. Run script 1 to initialize it.")
}

# --- PubChem Query Function for CID Properties ---
get_pubchem_properties <- function(cid) {
  url <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/", cid, "/property/SMILES,MolecularFormula,MonoisotopicMass/XML")
  
  tryCatch({
    xml <- read_xml(url)
    ns <- xml_ns(xml)
    
    list(
      SMILES = xml_text(xml_find_first(xml, ".//d1:SMILES", ns)),
      MolecularFormula = xml_text(xml_find_first(xml, ".//d1:MolecularFormula", ns)),
      MonoisotopicMass = as.numeric(xml_text(xml_find_first(xml, ".//d1:MonoisotopicMass", ns)))
    )
  }, error = function(e) {
    message(paste("  ERROR fetching properties for CID", cid, ":", e$message))
    return(NULL)
  })
}

# --- Enrich cache with missing values ---
for (i in seq_len(nrow(cid_cache_df))) {
  row <- cid_cache_df[i, ]
  
  if (!is.na(row$CID) && row$CID != 0 &&
      (is.na(row$SMILES) || is.na(row$MolecularFormula) || is.na(row$MonoisotopicMass))) {
    
    message(paste0("\n[", i, "/", nrow(cid_cache_df), "] Updating CID ", row$CID, " (", row$LookupName, ")"))
    props <- get_pubchem_properties(row$CID)
    Sys.sleep(0.2)
    
    if (!is.null(props)) {
      if ((is.na(row$SMILES) || row$SMILES == "") && !is.null(props$SMILES)) {
        message(paste("  ➤ SMILES updated:", props$SMILES))
        cid_cache_df$SMILES[i] <- props$SMILES
      }
      if ((is.na(row$MolecularFormula) || row$MolecularFormula == "") && !is.null(props$MolecularFormula)) {
        message(paste("  ➤ MolecularFormula updated:", props$MolecularFormula))
        cid_cache_df$MolecularFormula[i] <- props$MolecularFormula
      }
      if (is.na(row$MonoisotopicMass) && !is.null(props$MonoisotopicMass)) {
        message(paste("  ➤ MonoisotopicMass updated:", props$MonoisotopicMass))
        cid_cache_df$MonoisotopicMass[i] <- props$MonoisotopicMass
      }
      
    } else {
      message("  No data returned from PubChem.")
    }
  }
}

# --- Save updated cache ---
write.csv(cid_cache_df, cid_cache_file, row.names = FALSE)
message("\n✅ Cache updated and written to disk.")
