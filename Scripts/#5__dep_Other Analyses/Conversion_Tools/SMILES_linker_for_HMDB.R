##SMILES linker based on HMDB

#using json_df following from json_to_csv

standards <- read.csv("C18-Standard List Comparisons.csv")

json_df$HMDB_ID <- lapply(json_df$HMDB_ID, as.character)
standards$HMDB <- as.character(standards$HMDB)

library(dplyr)
library(tidyr)
library(stringr)

# Clean HMDB_ID values using regular expressions (keep as lists)
json_df <- json_df %>%
  mutate(HMDB_ID = lapply(HMDB_ID, function(x) str_extract_all(x, "HMDB\\d+"))) 

# Explode the HMDB_ID list column to have one row per HMDB ID per original row
json_df_long <- json_df %>% 
  unnest(HMDB_ID) %>%

  json_df_long <- json_df_long %>% 
rename(HMDB = 1) 

json_df_long$HMDB <- as.character( json_df_long$HMDB)

# Join with standards (using unnested json_df_long)
merged_df <- json_df_long %>%
  left_join(standards, by = "HMDB") 

# Group by the original row identifier (assuming you have one) and collapse smiles back into a list
result_df <- merged_df %>%
  group_by(row_number()) %>%  # Replace row_number() with your actual row identifier
  summarise(HMDB = first(HMDB),  # Rename HMDB_ID to HMDB
            smiles = list(smiles), .groups = "drop") 

# Join with the original standards df to add the smiles column
standards <- standards %>%
  left_join(result_df, by = "HMDB")

standards <- standards %>%
  mutate(smiles = sapply(smiles, paste, collapse = ",")) 
write.csv(standards, "standards.csv", row.names = FALSE)