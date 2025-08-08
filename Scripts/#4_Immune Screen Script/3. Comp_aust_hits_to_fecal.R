##Bioactivity correlator to immune cell screen
library(dplyr)

#Note: current deficiency in naming, need to get library codes from VCFG
CA.lib <- read.csv("~/CA/CompAus_Overlapping_CIDs.csv")

#results of overlap
A.unstim <- read.csv("~/PMC266/Microbiome McConville Lab Primary Screen PMC266  BMDM Subpopulation Percentages - Results   Plate A Unstimulated.csv")
A.unstim$Compound_Name <- gsub("_", " ", A.unstim$Compound_Name)
B.negpos <- read.csv("~/PMC266/Microbiome McConville Lab Primary Screen PMC266  BMDM Subpopulation Percentages - Results   Plate B IL4  IL13 Stimulated (1).csv")
B.negpos$Compound_Name <- gsub("_", " ", B.negpos$Compound_Name)
B.pospos <- read.csv("~/PMC266/Microbiome McConville Lab Primary Screen PMC266  BMDM Subpopulation Percentages - Results   Plate B IL4  IL13 Stimulated (2).csv")
B.pospos$Compound_Name <- gsub("_", " ", B.pospos$Compound_Name)
C.negpos <- read.csv("~/PMC266/Microbiome McConville Lab Primary Screen PMC266  BMDM Subpopulation Percentages - Results   Plate C LPS  IFNg Stimulated (1).csv")
C.negpos$Compound_Name <- gsub("_", " ", C.negpos$Compound_Name)

#total results
M0_pos_inos <- read.csv("~/PMC266/KW_PMC266_FC 15_Hit_Dose_Count_M0_pos_inos.csv")
M0_pos_inos$Compound_Name <- gsub("_", " ", M0_pos_inos$Compound_Name)
M0_pos_arg <- read.csv("~/PMC266/KW_PMC266_FC 15_Hit_Dose_Count_M0_pos_arg.csv")
M0_pos_arg$Compound_Name <- gsub("_", " ", M0_pos_arg$Compound_Name)

M1_pos_inos <- read.csv("~/PMC266/KW_PMC266_FC 0.2_Hit_Dose_Count_M1_pos_inos.csv")
M1_pos_inos$Compound_Name <- gsub("_", " ", M1_pos_inos$Compound_Name)
M1_pos_arg <- read.csv("~/PMC266/KW_PMC266_FC 15_Hit_Dose_Count_M1_pos_arg.csv")
M1_pos_arg$Compound_Name <- gsub("_", " ", M1_pos_arg$Compound_Name)

M2_pos_inos <- read.csv("~/PMC266/KW_PMC266_FC 15_Hit_Dose_Count_M2_pos_inos.csv")
M2_pos_inos$Compound_Name <- gsub("_", " ", M2_pos_inos$Compound_Name)
M2_pos_arg <- read.csv("~/PMC266/KW_PMC266_FC 1.15_Hit_Dose_Count_M2_pos_arg.csv")
M2_pos_arg$Compound_Name <- gsub("_", " ", M2_pos_arg$Compound_Name)


# --- 1. Add a 'Source_Screen' column to each screen dataframe ---
# This tags each compound with its origin before we combine them.
A.unstim$Source_Screen <- "A_Unstimulated"
B.negpos$Source_Screen <- "B_IL4_IL13_Stim_1"
B.pospos$Source_Screen <- "B_IL4_IL13_Stim_2"
C.negpos$Source_Screen <- "C_LPS_IFNg_Stim"

M0_pos_inos$Source_Screen <- "M0_Pos_inos"
M0_pos_arg$Source_Screen <- "M0_Pos_arg"
M1_pos_inos$Source_Screen <- "M1_Pos_inos"
M1_pos_arg$Source_Screen <- "M1_Pos_arg"
M2_pos_inos$Source_Screen <- "M2_Pos_inos"
M2_pos_arg$Source_Screen <- "M2_Pos_arg"

# --- 2. Combine the essential columns from all screens into one dataframe ---
# We only need the Compound_Name and the new Source_Screen column.
all_screens_long <- rbind(
  A.unstim[, c("Compound_Name", "Source_Screen")],
  B.negpos[, c("Compound_Name", "Source_Screen")],
  B.pospos[, c("Compound_Name", "Source_Screen")],
  C.negpos[, c("Compound_Name", "Source_Screen")],
  M0_pos_inos[, c("Compound_Name", "Source_Screen")],
  M0_pos_arg[, c("Compound_Name", "Source_Screen")],
  M1_pos_inos[, c("Compound_Name", "Source_Screen")],
  M1_pos_arg[, c("Compound_Name", "Source_Screen")],
  M2_pos_inos[, c("Compound_Name", "Source_Screen")],
  M2_pos_arg[, c("Compound_Name", "Source_Screen")]
)

# --- 3. Filter this combined list to find overlaps with CA.lib ---
# The %in% operator creates a TRUE/FALSE vector to keep only the rows
# where a compound name is also present in CA.lib.
overlapping_compounds_with_source <- all_screens_long[
  all_screens_long$Compound_Name %in% CA.lib$Compound_Name,
]

# --- 4. Clean up the results for display ---
# Remove any potential duplicates (e.g., if a compound appeared twice in the same file).
overlapping_compounds_with_source <- unique(overlapping_compounds_with_source)

# Order the results alphabetically by compound name for easier reading.
overlapping_compounds_with_source <- overlapping_compounds_with_source[
  order(overlapping_compounds_with_source$Compound_Name),
]

# Reset the row names for a clean output.
rownames(overlapping_compounds_with_source) <- NULL

# --- Display the detailed result ---
print("--- DETAILED VIEW: Overlapping compounds and the screen they were found in ---")
print(overlapping_compounds_with_source)


# --- (Optional) 5. Create a summarized view ---
# This creates a dataframe where each compound has only one row,
# and the 'Screens' column lists all its sources, collapsed into one string.
summarized_overlap <- aggregate(
  Source_Screen ~ Compound_Name,
  data = overlapping_compounds_with_source,
  FUN = function(x) paste(unique(x), collapse = ", ")
)

# --- Display the summarized result ---
print("--- SUMMARIZED VIEW: Each compound with all its source screens ---")
print(summarized_overlap)

write_csv(summarized_overlap, "~/PMC266/PMC_266_CompAust_Hits.csv")