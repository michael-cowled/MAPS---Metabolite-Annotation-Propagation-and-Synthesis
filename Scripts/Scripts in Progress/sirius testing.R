canopus.data <- paste0(folder, "/sirius-top100/canopus_structure_summary-100.tsv")
csi.data <- paste0(folder, "/sirius-top100/structure_identifications_top-100.tsv")
zodiac.data <- paste0(folder, "/sirius-top100/formula_identifications_top-100.tsv")


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
# Filter and keep only unique feature.ID (in the case where multiple annotations are provided)
csi.data <- csi.data %>%
  group_by(feature.ID) %>%
  filter(!(all(csi.confidence.score == 0))) %>%  # Remove groups where  all scores are 0
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

## 8. Calculation of Identification Probability (Wishart paper)
gnps.data <- gnps.data %>%
  mutate(gnps.ID.prob = 0)
gnps.data <- gnps.data %>%
  group_by(feature.ID) %>%
  mutate(
    n_above_thresh = sum(gnps.cosine.score >= 0.7),
    gnps.ID.prob = ifelse(gnps.cosine.score >= 0.7 & n_above_thresh > 0, 1 / n_above_thresh, 0)
  ) %>%
  ungroup() %>%
  select(-n_above_thresh)
gnps.data <- gnps.data %>%
  group_by(feature.ID) %>%
  arrange(desc(gnps.cosine.score)) %>%
  slice(1) %>%
  ungroup()  