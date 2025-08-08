library(readr)
library(dplyr)
library(ggplot2)
library(ggvenn)

##------------ Load and process data -------------##
# Load datasets
df1 <- read_csv("Y:/MA_BPA_Microbiome/Total-List-Of-Annotations.csv")
df2 <- read.delim("HMDB_Feces_Metabolites.TSV", header = TRUE, sep = "\t")
df2 <- df2 %>%
  rename(CID = pubchem_compound_id)

# Ensure CID columns are same type and name
df1$CID <- as.numeric(df1$CID)
df2$cid <- as.numeric(df2$CID)


# Combine all unique CIDs (including NA)
all_cids <- data.frame(CID = union(df1$CID, df2$CID)) %>%
  distinct()

# Add presence/absence info
all_cids <- all_cids %>%
  mutate(
    `Presence/Absence in fecal samples` = ifelse(!is.na(CID) & CID %in% df1$CID, "Y", "N"),
    `Presence/Absence in HMDB` = ifelse(!is.na(CID) & CID %in% df2$CID, "Y", "N")
  )

# Join back for inspection (optional)
merged_df_fecal <- left_join(all_cids, df1, by = "CID")
merged_df_HMDB <- left_join(all_cids, df2, by = "CID")
write_csv(merged_df_fecal, "merged_df_fecal.csv", row.names = FALSE, na = "")
write_csv(merged_df_HMDB, "merged_df_HMDB.csv", row.names = FALSE, na = "")

##------------ Count special cases -------------##
fecal_na_count <- sum(is.na(df1$CID))
fecal_zero_count <- sum(df1$CID == "0", na.rm = TRUE)
fecal_invalid <- fecal_na_count + fecal_zero_count

HMDB_na_count <- sum(is.na(df2$CID))
HMDB_zero_count <- sum(df2$CID == "0", na.rm = TRUE)
HMDB_invlaid <-HMDB_na_count + HMDB_zero_count

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
      CID %in% cid_df2 ~ "HMDB Only",
      TRUE ~ "Unknown"
    )
  )

# Set lists
fecal_cids <- overlap_df %>%
  filter(Source %in% c("Fecal Samples Only", "Both")) %>%
  pull(CID) %>% unique()

HMDB_cids <- overlap_df %>%
  filter(Source %in% c("HMDB Only", "Both")) %>%
  pull(CID) %>% unique()

cid_sets <- list(
  `Fecal Samples` = fecal_cids,
  `HMDB` = HMDB_cids
)

##------------ Plot Venn with annotations -------------##
p1 <- ggvenn(cid_sets,
             fill_color = c("#fc8d62", "#8da0cb"),
             stroke_size = 0.5,
             set_name_size = 4) +
  ggtitle("CID Overlap Between Fecal Samples and HMDB") +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", color = NA)) +
  annotate("text", x = -1.3, y = 0.6,
           label = paste0("Fecal annotation without CID: ", fecal_invalid),
           size = 3, hjust = 0) +
  annotate("text", x = 0.5, y = 0.6,
           label = paste0("HMDB without CID: ", HMDB_invlaid),
           size = 3, hjust = 0)

print(p1)
ggsave("venn_enriched.png", p1, width = 6, height = 4, dpi = 300)


###------------------matching with SMILES and names for those does not have valid CID---------------------###

# Get rows in HMDB with NA CID
df_HMDB_no_cid <- df2 %>% filter(is.na(CID))

# Get rows in Total annotations with NA or "0" CID
df_feces_no_cid <- df1 %>% filter(is.na(CID) | CID == "0")

# Optionally, save these to CSV
write_csv(df_HMDB_no_cid, "HMDB_No_CID.csv", row.names = FALSE, na = "")
write_csv(df_feces_no_cid, "Fecal_No_CID.csv", row.names = FALSE, na = "")


# Ensure relevant columns are character
df_feces_no_cid$smiles <- as.character(df_feces_no_cid$smiles)
df_HMDB_no_cid$SMILES <- as.character(df_HMDB_no_cid$SMILES)
df_feces_no_cid$compound.name <- as.character(df_feces_no_cid$compound.name)
df_HMDB_no_cid$name <- as.character(df_HMDB_no_cid$name)

# Step 1: Match by SMILES
matched_smiles <- intersect(df_feces_no_cid$smiles, df_HMDB_no_cid$SMILES)
matched_by_smiles <- matched_smiles[!is.na(matched_smiles) & matched_smiles != ""]

# Step 2: Match by annotation name (name vs compound.name)
matched_names <- intersect(df_feces_no_cid$compound.name, df_HMDB_no_cid$name)
matched_by_names <- matched_names[!is.na(matched_names) & matched_names != ""]

# Step 3: Combine
no_cid_match_name_smiles <- union(matched_by_smiles, matched_by_names)

# Step 4: Filter both datasets using the matched vector

# In fecal, filter if match in smiles or compound.name
no_cid_namesmiles_match_fecal <- df_feces_no_cid %>%
  filter(smiles %in% no_cid_match_name_smiles |
           compound.name %in% no_cid_match_name_smiles)

# In HMDB, filter if match in SMILES or name
no_cid_namesmiles_match_HMDB <- df_HMDB_no_cid %>%
  filter(SMILES %in% no_cid_match_name_smiles |
           name %in% no_cid_match_name_smiles)

# Optionally, save results
write_csv(no_cid_namesmiles_match_fecal, "no_cid_namesmiles_match_fecal.csv", row.names = FALSE, na = "")
write_csv(no_cid_namesmiles_match_HMDB, "no_cid_namesmiles_match_HMDB.csv", row.names = FALSE, na = "")


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
write_csv(df2_overlap, "HMDB_Overlapping_CIDs.csv", row.names = FALSE, na = "")

