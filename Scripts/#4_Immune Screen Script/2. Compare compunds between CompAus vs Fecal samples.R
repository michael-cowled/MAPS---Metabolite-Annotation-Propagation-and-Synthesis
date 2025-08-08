library(readr)
library(dplyr)
library(ggplot2)
library(ggvenn)

##------------ Load and process data -------------##
setwd("~/CA")
# Load datasets
df1 <- read_csv("Y:/MA_BPA_Microbiome/Total-List-Of-Annotations.csv")
df2 <- read_csv("~/CompAus_Final_with_Cleaned_Names.csv")
df2 <- df2 %>% 
  filter(grepl("DMSO", Solubility))

# Ensure CID columns are same type and name
df1$CID <- as.character(df1$CID)

# Combine all unique CIDs (including NA)
all_cids <- data.frame(CID = union(df1$CID, df2$CID)) %>%
  distinct()

# Add presence/absence info
all_cids <- all_cids %>%
  mutate(
    `Presence/Absence in fecal samples` = ifelse(!is.na(CID) & CID %in% df1$CID, "Y", "N"),
    `Presence/Absence in CompAus` = ifelse(!is.na(CID) & CID %in% df2$CID, "Y", "N")
  )

# Join back for inspection (optional)
merged_df_fecal <- left_join(all_cids, df1, by = "CID")
df2$CID <- as.character(df2$CID)
merged_df_CompAus <- left_join(all_cids, df2, by = "CID")
write_csv(merged_df_fecal, "merged_df_fecal.csv")
write_csv(merged_df_CompAus, "merged_df_CompAus.csv")

##------------ Count special cases -------------##
fecal_na_count <- sum(is.na(df1$CID))
fecal_zero_count <- sum(df1$CID == "0", na.rm = TRUE)
fecal_invalid <- fecal_na_count + fecal_zero_count

compaus_na_count <- sum(is.na(df2$CID))
compaus_zero_count <- sum(df2$CID == "0", na.rm = TRUE)
compaus_invlaid <-compaus_na_count + compaus_zero_count

##------------ Prepare for Venn -------------##
# Unique CID sets
cid_df1 <- unique(df1$CID)
cid_df2 <- unique(df2$CID)

# Categorize
overlap_df <- data.frame(
  CID = union(cid_df1, cid_df2)
) %>%
  mutate(
    Source = case_when(
      CID %in% cid_df1 & CID %in% cid_df2 ~ "Both",
      CID %in% cid_df1 ~ "Fecal Samples Only",
      CID %in% cid_df2 ~ "CompAus Only",
      TRUE ~ "Unknown"
    )
  )

# Set lists
fecal_cids <- overlap_df %>%
  filter(Source %in% c("Fecal Samples Only", "Both")) %>%
  pull(CID) %>% unique()

compaus_cids <- overlap_df %>%
  filter(Source %in% c("CompAus Only", "Both")) %>%
  pull(CID) %>% unique()

cid_sets <- list(
  `Fecal Samples` = fecal_cids,
  `CompAus` = compaus_cids
)

##------------ Plot Venn with annotations -------------##
p1 <- ggvenn(cid_sets,
             fill_color = c("#fc8d62", "#8da0cb"),
             stroke_size = 0.5,
             set_name_size = 4) +
  ggtitle("CID Overlap Between Fecal Samples and CompAus") +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", color = NA)) +
  annotate("text", x = -1.3, y = 0.6,
           label = paste0("Fecal annotation without CID: ", fecal_invalid),
           size = 3, hjust = 0) +
  annotate("text", x = 0.5, y = 0.6,
           label = paste0("CompAus without CID: ", compaus_invlaid),
           size = 3, hjust = 0)

print(p1)
ggsave("venn_enriched.png", p1, width = 6, height = 4, dpi = 300)


###------------------matching with SMILES and names for those does not have valid CID---------------------###

# Get rows in CompAus with NA CID
df_compaus_no_cid <- df2 %>% filter(is.na(CID))

# Get rows in Total annotations with NA or "0" CID
df_feces_no_cid <- df1 %>% filter(is.na(CID) | CID == "0")

# Optionally, save these to CSV
write_csv(df_compaus_no_cid, "CompAus_No_CID.csv")
write_csv(df_feces_no_cid, "Fecal_No_CID.csv")


# Ensure relevant columns are character
df_feces_no_cid$smiles <- as.character(df_feces_no_cid$smiles)
df_compaus_no_cid$Smiles <- as.character(df_compaus_no_cid$Smiles)
df_feces_no_cid$compound.name <- as.character(df_feces_no_cid$compound.name)
df_compaus_no_cid$Compound_Name <- as.character(df_compaus_no_cid$Compound_Name)

# Step 1: Match by SMILES
matched_smiles <- intersect(df_feces_no_cid$smiles, df_compaus_no_cid$Smiles)
matched_by_smiles <- matched_smiles[!is.na(matched_smiles) & matched_smiles != ""]

# Step 2: Match by annotation name (Compound_Name vs compound.name)
matched_names <- intersect(df_feces_no_cid$compound.name, df_compaus_no_cid$Compound_Name)
matched_by_names <- matched_names[!is.na(matched_names) & matched_names != ""]

# Step 3: Combine
no_cid_match_name_smiles <- union(matched_by_smiles, matched_by_names)

# Step 4: Filter both datasets using the matched vector

# In fecal, filter if match in smiles or compound.name
no_cid_namesmiles_match_fecal <- df_feces_no_cid %>%
  filter(smiles %in% no_cid_match_name_smiles |
           compound.name %in% no_cid_match_name_smiles)

# In CompAus, filter if match in SMILES or Compound_Name
no_cid_namesmiles_match_compaus <- df_compaus_no_cid %>%
  filter(Smiles %in% no_cid_match_name_smiles |
           Compound_Name %in% no_cid_match_name_smiles)

# Optionally, save results
write_csv(no_cid_namesmiles_match_fecal, "no_cid_namesmiles_match_fecal.csv")
write_csv(no_cid_namesmiles_match_compaus, "no_cid_namesmiles_match_compaus.csv")


##########---------------Ouptput overlapping CID----------------###
##------------ Extract and Save Overlapping CIDs -------------##

# Step 1: Get overlapping CIDs (present in both df1 and df2)
overlapping_cids <- intersect(cid_df1, cid_df2)

# Step 2: Filter df2 to only those rows with overlapping CIDs
df2_overlap <- df2 %>% filter(CID %in% overlapping_cids)
df2_overlap <- df2_overlap %>% filter(!is.na(CID))

repeats <- df2_overlap %>%
  group_by(CID) %>%
  filter(n() > 1)
df2_overlap <- df2_overlap %>%
  distinct(CID, .keep_all = TRUE)

# Step 3: Write filtered df2 to CSV
write_csv(df2_overlap, "CompAus_Overlapping_CIDs.csv")

