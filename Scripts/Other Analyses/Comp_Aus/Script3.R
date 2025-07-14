polar_list <- read.csv("polar_standards.csv")

# Check for duplicates in the key column 'cas' in merged_df
duplicated_cas_merged_df <- any(duplicated(merged_df$cas))

# Check for duplicates in the key column 'cas' in polar_list
duplicated_cas_polar_list <- any(duplicated(polar_list$cas))

# If duplicates exist in either data frame, remove them
if (duplicated_cas_merged_df) {
  merged_df <- merged_df[!duplicated(merged_df$cas), ]
}

if (duplicated_cas_polar_list) {
  polar_list <- polar_list[!duplicated(polar_list$cas), ]
}

# Now merge the two data frames
merged_df2 <- merge(merged_df, polar_list, by = "cas", all = TRUE)

# Assuming 'merged_df2' is your dataframe

# Create a new column with true/false indicating presence of "DMSO"
merged_df2$DMSO_present <- grepl("DMSO", merged_df2$Solubility)

#Filtering out DMSO insolubles:
merged_df2 <- filter(merged_df2, DMSO_present == TRUE)  ##6910 rather than 8459 compounds
#NOte: also that this step removes the non-DMSO soluble standards that might have correlated
#Can confirm that if you remove this filter you regain the 532 standards.
