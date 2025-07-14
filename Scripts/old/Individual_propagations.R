propagated.annotation.data <- full.annotation.data %>%
  mutate(Probably.Analogue.of = NA)  # Initialize new column with NA

# Filter the rows where either column contains the value 332 using dplyr
filtered_pairs <- gnps.cluster.pairs %>%
  filter(CLUSTERID1 == 134 | CLUSTERID2 == 134) %>%
  arrange(desc(Cosine))  # Sort by decreasing Cosine

# Extract paired values
paired_values <- filtered_pairs %>%
  mutate(paired_value = ifelse(CLUSTERID1 == 134, CLUSTERID2, CLUSTERID1)) %>%
  pull(paired_value)

# Function to get the first non-NA value from specified columns
get_first_non_na_value <- function(paired_value, summary.annotation.data) {
  result <- summary.annotation.data %>%
    filter(feature.ID == paired_value) %>%
    select(gnps.compound.name, csi.compound.name, ms2query.analogue.compound.name, canopus.NPC.pathway) %>%
    unlist() %>%
    na.omit() %>%
    .[1]
  if (is.null(result)) {
    result <- summary_data %>%
      filter(feature.ID == paired_value) %>%
      select(csi.smiles) %>%
      unlist() %>%
      na.omit() %>%
      .[1]
  }
  return(result)
}

# Find and update the value for feature.ID 332
result <- NA  # Initialize result
for (value in paired_values) {
  result <- get_first_non_na_value(value, full.annotation.data)
  if (!is.na(result)) {
    break  # Exit the loop if a valid result is found
  }
}

# Update the new column in propagated.annotation.data
propagated.annotation.data <- propagated.annotation.data %>%
  mutate(Probably.Analogue.of = ifelse(feature.ID == 134, result, Probably.Analogue.of))