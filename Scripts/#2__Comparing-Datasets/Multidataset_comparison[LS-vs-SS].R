library(dplyr)
library(tidyr)
library(readr)
library(ggplot2) # Ensure ggplot2 is loaded for plotting
library(ggvenn) # Only needed if you want the Venn diagram for the intra-dataset comparison

##-----------USER INPUT------------##

# List of dataset IDs for inter-dataset comparison
list.of.datasets.inter <- c("HGMD_0071+75_combined", "HGMD_0108")
metadata_file_HGMD_0125 <- "metadata.txt" # Assuming a separate metadata file for 0125

##-----------END USER INPUT------------##

# Function to tidy datasets (from Dataset_Compare(Inter).R)
dataset.tidier <- function(list.of.datasets) {
  dataset_list <- list()
  
  for (i in list.of.datasets) {
    print(i)
    file_path <- paste0("Y:/MA_BPA_Microbiome/Dataset-Annotations/", i, ".csv")
    
    if (file.exists(file_path)) {
      dataset <- read_csv(file_path) # Use read_csv for consistency and usually better performance
      dataset$dataset.id <- i
      dataset <- dataset %>%
        select(dataset.id, feature.usi, Samples, Best.Annotation, Best.Annotation.Smiles, 
               Best.Annotation.Confidence.Level, Best.Annotation.Compound.Class)  %>%
        distinct(Best.Annotation, .keep_all = TRUE)
      
      dataset$Best.Annotation.Confidence.Level <- as.numeric(as.character(dataset$Best.Annotation.Confidence.Level))
      
      dataset_list[[i]] <- dataset
    } else {
      warning(paste("File not found:", file_path))
    }
  }
  combined_datasets <- bind_rows(dataset_list)
  return(combined_datasets)
}

# 1. Process the first two datasets
datasets_inter <- dataset.tidier(list.of.datasets.inter)

# Rename dataset.id for clarity
datasets_inter <- datasets_inter %>%
  mutate(dataset.id = case_when(
    dataset.id == "HGMD_0071+75_combined" ~ "Phe-Hex-Large-Scale",
    dataset.id == "HGMD_0108" ~ "Phe-Hex-Small-Scale",
    TRUE ~ dataset.id # Keep other IDs as is if any
  ))

# 2. Process HGMD_0125 using metadata (adapted from Meta_Compare(Intra).R)
file_path_HGMD_0125 <- "Y:/MA_BPA_Microbiome/Dataset-Annotations/HGMD_0125.csv"
metadata_HGMD_0125 <- read_tsv(metadata_file_HGMD_0125)

# Read the HGMD_0125 annotations
HGMD_0125_annotations <- read_csv(file_path_HGMD_0125) %>%
  select(feature.usi, Samples, Best.Annotation, Best.Annotation.Smiles, 
         Best.Annotation.Confidence.Level, Best.Annotation.Compound.Class) %>%
  # Ensure confidence level is numeric
  mutate(Best.Annotation.Confidence.Level = as.numeric(as.character(Best.Annotation.Confidence.Level))) %>%
  distinct(Best.Annotation, .keep_all = TRUE)

# Filter metadata for groups
group_HILIC_Large_Scale <- filter(metadata_HGMD_0125, ATTRIBUTE_Group == "HILIC-Large-Scale")
group_HILIC_Small_Scale <- filter(metadata_HGMD_0125, ATTRIBUTE_Group == "HILIC-Small-Scale")

# Initialize new columns for presence in groups
HGMD_0125_annotations$ContainsHILIC_Large_Scale <- FALSE
HGMD_0125_annotations$ContainsHILIC_Small_Scale <- FALSE

# Populate the presence columns based on filenames
for (i in 1:nrow(group_HILIC_Large_Scale)) {
  filename <- group_HILIC_Large_Scale$filename[i]  
  HGMD_0125_annotations$ContainsHILIC_Large_Scale <- HGMD_0125_annotations$ContainsHILIC_Large_Scale | grepl(filename, HGMD_0125_annotations$Samples)
}

for (i in 1:nrow(group_HILIC_Small_Scale)) {
  filename <- group_HILIC_Small_Scale$filename[i]  
  HGMD_0125_annotations$ContainsHILIC_Small_Scale <- HGMD_0125_annotations$ContainsHILIC_Small_Scale | grepl(filename, HGMD_0125_annotations$Samples)
}

# Create tidy dataframes for HILIC groups
desired_cols <- c(
  "dataset.id",
  "feature.usi",
  "Samples",
  "Best.Annotation",
  "Best.Annotation.Smiles",
  "Best.Annotation.Confidence.Level",
  "Best.Annotation.Compound.Class"
)

df_HILIC_Large_Scale_tidy <- HGMD_0125_annotations %>%
  filter(ContainsHILIC_Large_Scale == TRUE) %>%
  mutate(dataset.id = "HILIC-Large-Scale") %>%
  select(all_of(desired_cols))

df_HILIC_Small_Scale_tidy <- HGMD_0125_annotations %>%
  filter(ContainsHILIC_Small_Scale == TRUE) %>%
  mutate(dataset.id = "HILIC-Small-Scale") %>%
  select(all_of(desired_cols))

# 3. Combine all four datasets
all_combined_datasets <- bind_rows(
  datasets_inter,
  df_HILIC_Large_Scale_tidy,
  df_HILIC_Small_Scale_tidy
)

# Ensure correct column names (already handled by select in dataset.tidier and manual selection)
if(length(colnames(all_combined_datasets)) >= 5){
  if(colnames(all_combined_datasets)[5] != "Best.Annotation.Smiles"){
    stop("The fifth column is expected to be named 'Best.Annotation.Smiles'. Please rename it accordingly.")
  }
} else {
  stop("The 'all_combined_datasets' dataframe does not have at least 5 columns")
}

# Get unique dataset IDs and confidence levels from the combined data
unique_datasets_all <- unique(all_combined_datasets$dataset.id)
unique_confidence_levels_all <- unique(all_combined_datasets$Best.Annotation.Confidence.Level)

# Initialize an empty list to store the results
result_list_final <- list()
result_counter_final <- 1

for (target_dataset in unique_datasets_all) {
  for (target_confidence in unique_confidence_levels_all) {
    
    target_data <- all_combined_datasets %>%
      filter(dataset.id == target_dataset, Best.Annotation.Confidence.Level == target_confidence)
    
    if (nrow(target_data) == 0) {
      next
    }
    
    target_data_distinct <- target_data %>%
      distinct(Best.Annotation.Smiles, .keep_all = TRUE)
    
    other_data <- all_combined_datasets %>%
      filter(dataset.id != target_dataset, Best.Annotation.Confidence.Level == target_confidence)
    
    if (nrow(other_data) > 0) {
      filtered_data <- target_data_distinct %>%
        anti_join(other_data, by = c("Best.Annotation.Smiles" = "Best.Annotation.Smiles"))
    } else {
      filtered_data <- target_data_distinct
    }
    
    result_list_final[[result_counter_final]] <- filtered_data
    result_counter_final <- result_counter_final + 1
  }
}

final_result_all <- bind_rows(result_list_final)

print(final_result_all)
write.csv(final_result_all, "filtered_datasets_all_four.csv", row.names = FALSE)

# Verification
dataset_summary_exclusive_all <- final_result_all %>%
  group_by(dataset.id, Best.Annotation.Confidence.Level) %>%
  summarize(unique_smiles_count = n(), .groups = "drop")
print(dataset_summary_exclusive_all)
write.csv(dataset_summary_exclusive_all, "unique_smiles_all_four.csv", row.names = FALSE)

# Count USIs per dataset.id and confidence.level
usi_counts_all <- all_combined_datasets %>%
  group_by(dataset.id, Best.Annotation.Confidence.Level) %>%
  summarize(total_usi = n_distinct(feature.usi), .groups = "drop")

write_csv(usi_counts_all, "usi_counts_all_four.csv")
print(usi_counts_all)

# Compound Class Plot (Unique Features)
filtered_datasets_plot <- final_result_all %>%
  filter(!(Best.Annotation.Compound.Class %in% c("None", "Others", "N/A", "NA", "")) & !is.na(Best.Annotation.Compound.Class))

plot_data_unique <- filtered_datasets_plot %>%
  group_by(Best.Annotation.Compound.Class, dataset.id) %>%
  summarise(Feature_Count = n(), .groups = "drop") %>%
  # Filter categories with less than 10 annotations
  filter(Feature_Count >= 3)

ggplot(plot_data_unique, aes(x = Feature_Count, y = Best.Annotation.Compound.Class, fill = dataset.id)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "Unique Feature Count by Superclass and Dataset (All Four, >= 10 Annotations)",
    x = "Number of Annotations",
    y = "Superclass",
    fill = "Dataset ID"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    legend.position = "right"
  ) -> p_unique

ggsave(
  filename = "Best_Annotation_Compound_Class_Annotations_counts_unique_all_four.svg",
  plot = p_unique,
  width = 10,
  height = 14,
  dpi = 300
)

# Compound Class Plot (Total Features) - from 'all_combined_datasets' before filtering for uniqueness
plot_data_total <- all_combined_datasets %>%
  filter(!(Best.Annotation.Compound.Class %in% c("None", "Others", "N/A", "NA", "")) & !is.na(Best.Annotation.Compound.Class)) %>%
  group_by(Best.Annotation.Compound.Class, dataset.id) %>%
  summarise(Feature_Count = n(), .groups = "drop") %>%
  # Filter categories with less than 10 annotations
  filter(Feature_Count >= 10)

ggplot(plot_data_total, aes(x = Feature_Count, y = Best.Annotation.Compound.Class, fill = dataset.id)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "Total Annotations by Superclass and Dataset (All Four, >= 3 Annotations)",
    x = "Number of Annotations",
    y = "Superclass",
    fill = "Dataset ID"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    legend.position = "right"
  ) -> p_total

ggsave(
  filename = "Best_Annotation_Compound_Class_Annotation_counts_total_all_four.svg",
  plot = p_total,
  width = 10,
  height = 14,
  dpi = 300
)
