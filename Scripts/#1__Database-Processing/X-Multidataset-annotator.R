#### ANNOTATION TABLE SCRIPT FOR AUSTRALIAN HUMAN GUT METABOLOME DATABASE ####

#!# Run 1-Functions first! #!!#

#------------------------------------------------------------------------------#
###---START OF USER-FED INFORMATION---###

# Dataset IDs: From Data Management Plan:
# Define your list of dataset IDs here
dataset_ids_to_process <- good_datasets # Add all your desired IDs

# Specify the path to the Data Management Plan (unless moved to Mediaflux then should be consistent)
excel_file <- "C:/Users/mcowled/The University of Melbourne/MA Human Gut Metabolome - Documents/Data Management and Analysis/HGM - Data Management System.xlsx"

# Annotation Acceptance Probabilities (Defaults for Exploratory Analyses)
# Increase stringency if the application demands it.
gnps.prob <- 0.7
canopus.prob <- 0.7
csi.prob <- 0.64
ms2query.prob <- 0.63

###---END OF USER-FED INFORMATION---###
#------------------------------------------------------------------------------#

## 1. check_and_install
setwd(file.path("~/"))
# Function to check, install, and load required packages
check_and_install <- function(packages, github_packages = list()) {
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
                       "openxlsx", "tidyverse", "rvest", "jsonlite", "xml2", "progress")

github_packages <- list("MAPS.Package" = "michael-cowled/MAPS-Package-Public")
check_and_install(required_packages, github_packages)

# --- Load or Initialize Cache ---
# The cache should ideally be loaded once before the loop, and saved after the loop
# if you want to reuse it across all dataset processing.
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

# --- START LOOP FOR EACH DATASET ID ---
for (current_dataset_id in dataset_ids_to_process) {
  dataset.id <- current_dataset_id # Assign the current ID to the variable used in the script
  message(paste0("Processing Dataset ID: ", dataset.id))
  
  #-----------------------------------------------------------------------------------------------------------------------#
  ## 2. Updating Metadata and extract gnps task ID
  # This part reads from the excel file and creates CSVs.
  # Consider if this needs to be run for *every* dataset ID or only once if the Excel file is static.
  # If it's static and only dataset.id filter changes, move the `for (sheet in sheet_names)` loop outside.
  # If the Excel content can change per dataset process, keep it inside.
  # For now, I'll assume the CSVs are static or updated based on `excel_file` once.
  # For simplicity, I'll keep it in the loop but you might optimize this.
  sheet_names <- excel_sheets(excel_file)
  for (sheet in sheet_names) {
    data <- read_excel(excel_file, sheet = sheet)
    write.csv(data, file = paste0("HGM/", sheet, ".csv"), row.names = FALSE)
  }
  
  dataset <- read.csv("HGM/D - Dataset.csv") %>%
    filter(HGMD.ID == dataset.id)
  gnps.task.id <- dataset$gnps.task.ID[1]
  
  if (!is.na(dataset$column.type[1]) && dataset$column.type[1] == "HILIC") {
    rt.tol <- 0.2
  } else {
    rt.tol <- 0.1
  }
  
  if (is.na(gnps.task.id)) {
    stop(paste("gnps.task.ID is missing for dataset:", dataset.id, ". Add to 'dataset' csv before re-running."))
  }
  #-----------------------------------------------------------------------------------------------------------------------#
  
  ## 3. File List Extractor
  # ... (Rest of the code for step 3) ...
  dataset <- read.csv("HGM/D - Dataset.csv") %>%
    filter(HGMD.ID == dataset.id)
  output_directory <- dataset$Processed.Data.Folder[1]
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
    output_path <- file.path(output_directory, output_filename)
    write.csv(file_info, file = output_path, row.names = FALSE)
    print(paste("File list saved to:", output_path))
  } else {
    print("Folder does not exist. Check the path.")
  }
  
  ## 4. Processed Data Check
  # ... (Rest of the code for step 4, assuming `validate_and_get_paths` is defined in "1-Functions") ...
  folder <- dataset$Processed.Data.Folder[1]
  paths <- validate_and_get_paths(folder)
  mzmine.data <- paths$mzmine_data
  mzmine.annotations <- paths$mzmine_annotations
  canopus.data <- paths$canopus_data
  csi.data <- paths$csi_data
  zodiac.data <- paths$zodiac_data
  ms2query.data <- paths$ms2query_data
  cytoscape_file_path <- paths$cytoscape # Renamed to avoid conflict if `cytoscape` is a dataframe later
  
  ## 5. Load in MZMINE data
  # ... (Rest of the code for step 5, ensure `standardise_annotation` is defined in "1-Functions") ...
  mzmine.annotations <- read.csv(mzmine.annotations) %>%
    distinct(id, compound_name, .keep_all = TRUE) %>%
    group_by(id) %>%
    filter(score > min(score) | n() == 1) %>%
    ungroup() %>%
    group_by(compound_name) %>%
    filter(score == max(score)) %>%
    slice(1) %>%
    ungroup()
  
  mzmine.annotations$mzmine.id.prob <- NA
  mzmine.annotations <- standardise_annotation(mzmine.annotations, "compound_name", "smiles")
  # `cid_cache_df` is now global, so updates will persist across datasets if written after the loop.
  # write.csv(cid_cache_df, cid_cache_file, row.names = FALSE) # Move this outside the loop if you want one cumulative cache
  
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
    mutate(mzmine.id.prob = 1 / n()) %>%
    filter(score == max(score)) %>%
    slice(1) %>%
    ungroup() %>%
    select(id, compound_name, score, smiles, mzmine.id.prob, mol_formula, CID, MonoisotopicMass)
  names(mzmine.annotations.final) <- c('feature.ID', "compound.name", "confidence.score",
                                       "smiles", "id.prob", "molecular.formula", "CID", "MonoisotopicMass")
  mzmine.annotations.final$feature.ID <- as.numeric(mzmine.annotations.final$feature.ID)
  mzmine.annotations.final$confidence.level <- "1"
  mzmine.annotations.final$annotation.type <- "authentic standard"
  
  mzmine.data <- read.csv(mzmine.data)
  sample.data <- mzmine.data
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
  
  
  ## 6. Load, tidy and standardise GNPS2 data
  # ... (Rest of the code for step 6, ensure `fix_compound_names` is defined in "1-Functions") ...
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
  gnps.annotation.data <- gnps.annotation.data[, c(1, 6, 7, 3:5, 2, 8, 10, 9, 11, 12)]
  gnps.annotation.data <- gnps.annotation.data %>%
    filter(mz.diff.ppm <= 5) %>%
    filter(confidence.score >= gnps.prob)
  
  gnps.cluster.data <- read_tsv(paste0("https://gnps2.org/resultfile?task=", gnps.task.id, "&file=nf_output/networking/clustersummary_with_network.tsv"))
  gnps.cluster.data <- select(gnps.cluster.data, 'cluster index', 'component')
  names(gnps.cluster.data) <- c('feature.ID', "gnps.cluster.ID")
  gnps.cluster.pairs <- read_tsv(paste0("https://gnps2.org/resultfile?task=", gnps.task.id, "&file=nf_output/networking/filtered_pairs.tsv"))
  gnps.annotation.data$annotation.type <- "gnps"
  
  unique_in_gnps <- setdiff(gnps.annotation.data$feature.ID, mzmine.annotations.final$feature.ID)
  gnps.annotation.data <- filter(gnps.annotation.data, feature.ID %in% unique_in_gnps)
  
  gnps.data.lv2 <- filter(gnps.annotation.data, library.quality != "Insilico")
  gnps.data.lv2$confidence.level <- "2"
  gnps.data.lv2 <- fix_compound_names(gnps.data.lv2, "compound.name")
  
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
  # ... (Rest of the code for step 7) ...
  ms2query.data <- read.csv(ms2query.data)
  ms2query.data  <- ms2query.data[, c(2, 3, 4, 7, 8, 10, 17, 18)]
  names(ms2query.data) <- c("confidence.score", "mz.diff", "precursor_mz",
                            "compound.name", "smiles",
                            'feature.ID', "NPC.superclass", "NPC.pathway")
  ms2query.data <- fix_compound_names(ms2query.data, "compound.name")
  ms2query.data <- mutate(ms2query.data, mz.diff.ppm = mz.diff/precursor_mz *1000000) %>%
    filter(mz.diff.ppm <= 5)
  ms2query.data$annotation.type <- "ms2query"
  
  unique_in_ms2query <- setdiff(ms2query.data$feature.ID, mzmine.annotations.final$feature.ID)
  ms2query.data <- filter(ms2query.data, feature.ID %in% unique_in_ms2query)
  
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
  
  ms2query.data.lv3 <- ms2query.data %>%
    filter(confidence.score > ms2query.prob) %>%
    filter(mz.diff >= 0.001) %>%
    filter(!feature.ID %in% mzmine.annotations.final$feature.ID)
  ms2query.data.lv3$confidence.level <- "3"
  
  
  ## 8. Standardise level 2 annotations and compute ID probability
  # ... (Rest of the code for step 8, ensure `deduplicate_data` is defined in "1-Functions") ...
  lv2.annotations <- rbind(gnps.data.lv2, ms2query.data.lv2)
  
  lv2.annotations$CID <- NA
  lv2.annotations$MolecularFormula<- NA
  lv2.annotations$MonoisotopicMass <- NA
  lv2.annotations <- deduplicate_data(lv2.annotations, compound.name, confidence.score)
  lv2.annotations<- standardise_annotation(lv2.annotations, "compound.name", "smiles")
  lv2.annotations <- deduplicate_data(lv2.annotations, compound.name, confidence.score)
  # write.csv(cid_cache_df, cid_cache_file, row.names = FALSE) # Move this outside the loop if you want one cumulative cache
  
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
    ungroup()
  
  missing_cols <- colnames(lv2.annotations)[!colnames(lv2.annotations) %in% colnames(mzmine.annotations.final)]
  for (col in missing_cols) {
    mzmine.annotations.final[[col]] <- NA
  }
  mzmine.annotations.final <- select(mzmine.annotations.final, -molecular.formula)
  
  lv1.and.lv2.annotations <- lv2.annotations %>%
    rbind(mzmine.annotations.final)
  
  
  ## 9. Appending lv3 In silico matches from GNPS
  # ... (Rest of the code for step 9) ...
  unique_in_lv3<- setdiff(gnps.data.lv3$feature.ID, lv1.and.lv2.annotations$feature.ID)
  gnps.data.lv3 <- filter(gnps.data.lv3, feature.ID %in% unique_in_lv3)
  
  lv1.lv2.lv3.annotations <- lv1.and.lv2.annotations %>%
    rbind(gnps.data.lv3)
  
  ## 10. Load in SIRIUS data
  # ... (Rest of the code for step 10) ...
  canopus.data <- read_tsv(canopus.data)
  canopus.data <- canopus.data[, c(5:8, 27)]
  names(canopus.data) <- c("canopus.NPC.pathway", "canopus.NPC.pathway.probability",
                           "canopus.NPC.superclass", "canopus.NPC.superclass.probability",
                           'feature.ID')
  canopus.data <- canopus.data %>%
    group_by(feature.ID) %>%
    filter(!(all(canopus.NPC.pathway.probability == 0))) %>%
    filter(canopus.NPC.pathway.probability == max(canopus.NPC.pathway.probability, na.rm = TRUE)) %>%
    slice(1) %>%
    ungroup()
  
  zodiac.data <- read_tsv(zodiac.data)
  zodiac.data  <- zodiac.data[, c(2, 5, 21)]
  names(zodiac.data) <- c("zodiac.formula", "zodiac.confidence.score", 'feature.ID')
  zodiac.data <- zodiac.data %>%
    group_by(feature.ID) %>%
    filter(!(all(zodiac.confidence.score == 0))) %>%
    filter(zodiac.confidence.score == max(zodiac.confidence.score, na.rm = TRUE)) %>%
    slice(1) %>%
    ungroup()
  
  csi.data <- read_tsv(csi.data)
  csi.data <- csi.data[, c(3, 14, 15, 25)]
  names(csi.data) <- c("confidence.score", "compound.name", "smiles", 'feature.ID')
  csi.data <- csi.data[, c(2, 1, 3:ncol(csi.data))]
  csi.data$confidence.score[csi.data$confidence.score == -Inf] <- 0
  csi.data$confidence.score <- as.numeric(csi.data$confidence.score)
  csi.data <- filter(csi.data, confidence.score >= csi.prob)
  csi.data$compound.name[grepl("Solaparnaine", csi.data$compound.name, ignore.case = TRUE)] <- "Solaparnaine"
  csi.data$compound.name[grepl("Spectalinine", csi.data$compound.name, ignore.case = TRUE)] <- "(-)-Spectalinine"
  
  
  ## 11. Appending lv3 In silico matches from CSI:FingerID
  # ... (Rest of the code for step 11) ...
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
  
  unique_in_csi <- setdiff(csi.data$feature.ID, lv1.lv2.lv3.annotations$feature.ID)
  csi.data <- filter(csi.data, feature.ID %in% unique_in_csi)
  
  missing_cols <- colnames(lv1.lv2.lv3.annotations)[!colnames(lv1.lv2.lv3.annotations) %in% colnames(csi.data)]
  for (col in missing_cols) {
    csi.data[[col]] <- NA
  }
  csi.data$annotation.type <- "CSI:FingerID"
  csi.data$confidence.level <- "3"
  
  lv1.lv2.lv3.annotations <- lv1.lv2.lv3.annotations %>%
    rbind(csi.data)
  
  ## 12. Appending lv3 analogues from ms2query
  # ... (Rest of the code for step 12) ...
  ms2query.data.lv3$compound.name <- paste0("Analogue of ", ms2query.data.lv3$compound.name)
  
  unique_in_ms2query <- setdiff(ms2query.data.lv3$feature.ID, lv1.lv2.lv3.annotations$feature.ID)
  ms2query.data.lv3 <- filter(ms2query.data.lv3, feature.ID %in% unique_in_ms2query) %>%
    select(-mz.diff, -precursor_mz)
  
  missing_cols <- colnames(lv1.lv2.lv3.annotations)[!colnames(lv1.lv2.lv3.annotations) %in% colnames(ms2query.data.lv3)]
  for (col in missing_cols) {
    ms2query.data.lv3[[col]] <- NA
  }
  
  lv1.lv2.lv3.annotations <- lv1.lv2.lv3.annotations %>%
    rbind(ms2query.data.lv3)
  
  ## 13. Appending all other features
  # ... (Rest of the code for step 13) ...
  full.annotation.data <- mzmine.data[,1:4] %>%
    full_join(lv1.lv2.lv3.annotations, by = "feature.ID") %>%
    full_join(canopus.data, by = "feature.ID") %>%
    full_join(zodiac.data, by = "feature.ID") %>%
    full_join(gnps.cluster.data, by = "feature.ID")
  
  full.annotation.data$feature.usi <- paste0("mzspec:GNPS2:TASK-", gnps.task.id, "-nf_output/clustering/spectra_reformatted.mgf:scan:", full.annotation.data$feature.ID)
  
  ## 14. Propagation of annotations
  # ... (Rest of the code for step 14, ensure `paired_feature_finder` and `get_result` are defined in "1-Functions") ...
  na.rows <- filter(full.annotation.data, is.na(compound.name))
  na.feature.ids <- na.rows$feature.ID
  propagated_results <- list()
  
  pb <- progress::progress_bar$new(
    format = "Propagating annotations [:bar] :percent eta: :eta",
    total = length(na.feature.ids),
    width = 60
  )
  
  for (i in na.feature.ids) {
    pb$tick()
    paired_values <- paired_feature_finder(i)
    selected_paired_value <- NA
    final_result_data <- list(value = NA, column = NA, superclass = NA)
    
    for (value in paired_values) {
      result_data <- get_result(value, full.annotation.data)
      if (!is.na(result_data$value)) {
        selected_paired_value <- value
        final_result_data <- result_data
        annotation_type_value <- full.annotation.data %>%
          filter(feature.ID == selected_paired_value) %>%
          pull(annotation.type) %>%
          first()
        final_result_data$column <- annotation_type_value
        break
      }
    }
    
    if (!is.na(selected_paired_value)) {
      propagated_results[[as.character(i)]] <- list(
        feature.ID = i,
        Probable.Analogue.Of = final_result_data$value,
        Propagated.Feature.ID = selected_paired_value,
        Propagated.Annotation.Type = final_result_data$column,
        Propagated.Annotation.Class = final_result_data$superclass
      )
    }
  }
  
  propagated_df <- bind_rows(propagated_results)
  
  propagated.annotation.data <- full.annotation.data %>%
    left_join(propagated_df, by = "feature.ID")
  
  ## 15. Append propagations
  # ... (Rest of the code for step 15) ...
  propagated.annotation.data <- propagated.annotation.data %>%
    mutate(
      confidence.level = as.character(confidence.level),
      propagation_mask = is.na(compound.name) & !is.na(Probable.Analogue.Of),
      compound.name = ifelse(
        propagation_mask,
        paste0("Probable analogue of: ", Probable.Analogue.Of),
        compound.name
      ),
      confidence.level = ifelse(
        propagation_mask,
        "3",
        confidence.level
      ),
      NPC.pathway = NA,
      NPC.superclass = Propagated.Annotation.Class,
      annotation.type = ifelse(
        propagation_mask,
        "GNPS Propagation",
        annotation.type
      )
    ) %>%
    select(-propagation_mask)
  
  ## 16. Append level 4 and 5 annotations
  # ... (Rest of the code for step 16) ...
  propagated.annotation.data <- propagated.annotation.data %>%
    mutate(
      confidence.level = as.character(confidence.level),
      level4_mask = is.na(compound.name) &
        !is.na(canopus.NPC.pathway) &
        canopus.NPC.pathway.probability >= canopus.prob,
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
  
  
  ## 17. Sample List Data
  # ... (Rest of the code for step 17) ...
  sample.data <- sample.data %>%
    select(id, contains(".area"))
  
  colnames(sample.data)[1] <- "feature.ID"
  colnames_sample <- colnames(sample.data)
  colnames_sample <- sub("^datafile\\.", "", colnames_sample)
  colnames_sample <- sub("\\.mzML\\.area$", "", colnames_sample)
  colnames(sample.data) <- colnames_sample
  
  sample.data2 <- sample.data
  
  sample.data[, 2:ncol(sample.data)] <- lapply(sample.data[, 2:ncol(sample.data)], function(x) {
    x[x > 0] <- 1
    return(x)
  })
  
  sample.data$Samples <- apply(sample.data[, 2:ncol(sample.data)], 1, function(row) {
    paste(colnames(sample.data[, 2:ncol(sample.data)])[which(row == 1)], collapse = "; ")
  })
  
  sample.data <- sample.data[, c(1, ncol(sample.data))]
  
  propagated.annotation.data.with.samples <- propagated.annotation.data %>%
    full_join(sample.data, by = "feature.ID")
  
  ## 18. Cytoscape Appending
  # ... (Rest of the code for step 18) ...
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
  
  
  ## 19. Collapsing Ion Identity Networks
  # ... (Rest of the code for step 19, ensure `process_all_features` and `update_data_frame` are defined in "1-Functions") ...
  final.annotation.df <- propagated.annotation.data.with.samples
  new.summary.annotation.data <- propagated.annotation.data.with.samples
  iin.features <- filter(new.summary.annotation.data, !is.na(ion.identity.ID))
  
  # The `data` object for `process_all_features` seems to be undefined in the provided snippet here.
  # Assuming 'data' refers to `new.summary.annotation.data` or a specific subset of it.
  # For safety, I'll assume it's `new.summary.annotation.data` based on context.
  results <- process_all_features(new.summary.annotation.data, "ion.identity.ID", c("feature1", "feature2"))
  
  for (id in names(results)) {
    cat("Results for ion.identity.ID =", id, ":\n")
    print(results[[id]])
  }
  
  final.annotation.df <- update_data_frame(final.annotation.df, results, "ion.identity.ID") %>%
    select(-ion.identity.ID, -Probable.Analogue.Of)
  
  if ("Samples" %in% names(sample.data2)) {
    sample.data2 <- sample.data2 %>% select(-Samples)
  }
  
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
  
  final.annotation.df <- fix_compound_names(final.annotation.df, "compound.name")
  
  
  ## 20. Redundancy reduction
  # ... (Rest of the code for step 20, ensure `redundancy_fixer` is defined in "1-Functions") ...
  dataset_for_redundancy <- final.annotation.df %>%
    filter(!is.na(smiles) & smiles != "N/A")
  dataset_for_redundancy$rt <- as.numeric(dataset_for_redundancy$rt)
  dataset_for_redundancy$redundant <- FALSE
  
  other.data <- final.annotation.df %>%
    filter(is.na(smiles) | smiles == "N/A")
  
  dataset_for_redundancy <- redundancy_fixer(dataset_for_redundancy, column_to_check = "smiles")
  
  ## 21. Writing files to disk
  ################################
  ##ANNOTATION FILE
  final.annotation.df2 <- dataset_for_redundancy %>%
    filter(!redundant) %>%
    select(-redundant) %>%
    rbind(other.data)
  
  # Adjust output paths to include dataset.id for unique files
  write.csv(final.annotation.df2, paste0(folder, "/", dataset.id, "-final-annotation-df.csv"), row.names = FALSE)
  write.csv(final.annotation.df2, paste0("Y:/MA_BPA_Microbiome/Dataset-Annotations/", dataset.id, ".csv"), row.names = FALSE)
  
  ##################################
  ##SAMPLES DF
  output_file <- paste0(folder, "/", dataset.id, "-samples-df.csv") # Include folder and dataset.id
  chunk_size <- 50000
  total_rows <- nrow(samples.df)
  num_chunks <- ceiling(total_rows / chunk_size)
  
  write.table(samples.df[0, ], file = output_file, sep = ",",
              col.names = TRUE, row.names = FALSE, append = FALSE,
              qmethod = "double")
  
  cat(paste0("Starting to write to CSV for ", dataset.id, "...\n"))
  
  progress_bar <- txtProgressBar(min = 0, max = num_chunks, style = 3)
  
  for (i in 1:num_chunks) {
    start_row <- (i - 1) * chunk_size + 1
    end_row <- min(i * chunk_size, total_rows)
    
    chunk_data <- samples.df[start_row:end_row, ]
    
    write.table(chunk_data, file = output_file, sep = ",",
                col.names = FALSE, row.names = FALSE, append = TRUE,
                qmethod = "double")
    
    setTxtProgressBar(progress_bar, i)
  }
  close(progress_bar)
  cat(sprintf("\nCSV file saved to: %s\n", output_file))
  
  final.annotation.df3 <- final.annotation.df2 %>%
    select(feature.ID, mz, compound.name, smiles, confidence.level)
  
  samples.df.with.annotations <- samples.df %>%
    full_join(final.annotation.df3, by = "feature.ID")
  
  top_10_features <- samples.df.with.annotations %>%
    group_by(samples) %>%
    top_n(10, area) %>%
    ungroup()
  
  write.csv(top_10_features, paste0(folder, "/", dataset.id, "-top-10-features.csv"), row.names = FALSE) # Include folder and dataset.id
  
  message(paste0("Script finished processing for Dataset ID: ", dataset.id, "\n"))
} # --- END LOOP FOR EACH DATASET ID ---

# Save the updated CID cache after all datasets are processed
write.csv(cid_cache_df, cid_cache_file, row.names = FALSE)

print("All specified datasets have been processed.")