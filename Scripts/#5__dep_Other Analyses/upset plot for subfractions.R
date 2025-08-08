# Load the necessary libraries
library(tidyverse)
library(UpSetR)

# Load the dataset
metabolomics_data<- read.csv("Y:/MA_BPA_Microbiome/Dataset-Annotations/HGMD_0096.csv")

# --- Data Preparation ---
# 1. Select the feature ID and the Samples columns.
# 2. Split the semicolon-separated 'Samples' string into a list of individual sample names.
# 3. Transform the data into a "long" format, where each row represents one feature in one sample.
# 4. Group by sample name and create a list of all feature IDs for each sample.
# 5. Convert the result into a named list, which is the required format for UpSetR.

input_list <- metabolomics_data %>%
  select(feature.ID, Samples) %>%
  mutate(SampleList = str_split(Samples, ";\\s*")) %>%
  select(feature.ID, SampleList) %>%
  unnest(SampleList) %>%
  group_by(SampleList) %>%
  summarise(features = list(feature.ID), .groups = 'drop') %>%
  deframe()

# --- Generate the UpSet Plot ---
# Create the plot. You can adjust 'nsets' to show more or fewer samples.
# 'order.by = "freq"' sorts the intersections by size (most frequent first).
upset(
  fromList(input_list),
  nsets = 20, # Show intersections for the 20 most frequent sets (samples)
  order.by = "freq",
  text.scale = 1.5 # Adjust text size for better readability
)