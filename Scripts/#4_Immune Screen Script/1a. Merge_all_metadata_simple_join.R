###---------------------------Bind all data together---------------------------------------######
## Load data and clean column names  ###
#### New folder has been created called HGM_data_merge_metadata  #####
setwd("~/HGM")

library(dplyr)
library(purrr)
library(tidyr)
library(readxl)
library(tibble)

# Load CSVs
fractions_df <- read.csv("F - Fractions.csv", stringsAsFactors = FALSE)
analysis_df <- read.csv("A - Analysis.csv", stringsAsFactors = FALSE) %>%
  mutate(Injection.Volume.uL = as.character(Injection.Volume.uL))
extracts_df <- read.csv("E - Extracts.csv", stringsAsFactors = FALSE)
dataset_df <- read.csv("D - Dataset.csv", stringsAsFactors = FALSE)
     dataset_df <- dataset_df[, -c(15, 16, 17, 18, 19)]
immune_df <- read.csv("I - Immune Screen Samples.csv", stringsAsFactors = FALSE)



### Step 1: Construct column.info strings ###
fractions_df <- fractions_df %>%
  left_join(column_info_df, by = c("Column.ID" = "HGMC.ID")) 

dataset_df <- dataset_df %>%
  left_join(column_info_df, by = c("column.ID" = "HGMC.ID")) 


### Step 3: Prefix columns ### Exceptfor the Fractions df because the prefix was added after subsetting  ####
prefix_cols <- function(df, prefix) {
  names(df) <- paste0(prefix, names(df))
  df
}

analysis_df <- prefix_cols(analysis_df, "A.")
extracts_df <- prefix_cols(extracts_df, "E.")
fractions_df<-prefix_cols(fractions_df, "F.")
dataset_df <- prefix_cols(dataset_df, "D.")
immune_df <- prefix_cols(immune_df, "I.")

### Step 4: Determine parent type ###
analysis_df <- analysis_df %>%
  mutate(A.Parent.Type = case_when(
    grepl("^HGME", A.Parent.Code) ~ "Extract",
    grepl("^HGMF", A.Parent.Code) ~ "Fraction",
    TRUE ~ NA_character_
  ))




### Step 5: Join metadata ###
# Step 1: Base join with Analysis
joined_df <- analysis_df %>%
  left_join(dataset_df, by = c("A.dataset.ID" = "D.HGMD.ID")) %>%
  # F0: Direct parent is a fraction (non-fraction parent)
  left_join(fractions_df, by = c("A.Parent.Code" = "F.HGMF.ID")) %>%
  left_join(immune_df, by = c("A.Parent.Code" = "I.Fraction.Code"))



### Step 7: Export ###
write_csv(joined_df, "full_metadata_prefixed_joined_simple_join.csv", row.names = FALSE)






