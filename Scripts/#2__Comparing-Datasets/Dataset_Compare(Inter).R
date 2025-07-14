library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

##-----------USER INPUT------------##

# List of dataset IDs
list.of.datasets <- c("HGMD_0105", "HGMD_0108") # Add as many or little as you want.

##-----------END USER INPUT------------##

dataset.tidier <- function(list.of.datasets) {
  dataset_list <- list()
  
  for (i in list.of.datasets) {
    print(i)
    file_path <- paste0("Y:/MA_BPA_Microbiome/Dataset-Annotations/", i, ".csv")
    
    if (file.exists(file_path)) {
      dataset <- read.csv(file_path)
      dataset$dataset.id <- i
      dataset <- dataset %>%
        select(dataset.id, feature.usi, Samples, Best.Annotation, Best.Annotation.Smiles, 
               Best.Annotation.Confidence.Level, Best.Annotation.Compound.Class)      # %>%
      # Optional filtering step (remove # above and below)
        # distinct(Best.Annotation, .keep_all = TRUE)
      
      # --- FIX START ---
      # Convert Best.Annotation.Confidence.Level to numeric, coercing non-numeric to NA
      dataset$Best.Annotation.Confidence.Level <- as.numeric(as.character(dataset$Best.Annotation.Confidence.Level))
      # --- FIX END ---
      
      dataset_list[[i]] <- dataset
    } else {
      warning(paste("File not found:", file_path))
    }
  }
  combined_datasets <- bind_rows(dataset_list)
  return(combined_datasets)
}

# Call the function and store the result
datasets <- dataset.tidier(list.of.datasets)
  
# Ensure correct column names
if(length(colnames(datasets))>=5){
  if(colnames(datasets)[5] != "Best.Annotation.Smiles"){
    stop("The fifth column is expected to be named 'Best.Annotation.Smiles'. Please rename it accordingly.")
  }
} else {
  stop("The 'datasets' dataframe does not have at least 5 columns")
}


# Get unique dataset IDs and confidence levels
unique_datasets <- unique(datasets$dataset.id)
unique_confidence_levels <- unique(datasets$Best.Annotation.Confidence.Level)

# Initialize an empty list to store the results
result_list <- list()
result_counter <- 1  # To keep track of the list index

for (target_dataset in unique_datasets) {
  for (target_confidence in unique_confidence_levels) {
    
    # 1. Subset the data for the target dataset and confidence level
    target_data <- datasets %>%
      filter(dataset.id == target_dataset, Best.Annotation.Confidence.Level == target_confidence)
    
    # Skip to the next iteration if target_data is empty
    if (nrow(target_data) == 0) {
      next
    }
    
    # 2. Remove duplicates *within* the target dataset/confidence
    target_data_distinct <- target_data %>%
      distinct(Best.Annotation.Smiles, .keep_all = TRUE)
    
    # 3. Subset the data for all *other* datasets and the *same* confidence level
    other_data <- datasets %>%
      filter(dataset.id != target_dataset, Best.Annotation.Confidence.Level == target_confidence)
    
    # 4. Perform the comparison and filtering (efficiently with anti_join)
    if (nrow(other_data) > 0) {
      filtered_data <- target_data_distinct %>%
        anti_join(other_data, by = c("Best.Annotation.Smiles" = "Best.Annotation.Smiles"))
    } else {
      filtered_data <- target_data_distinct
    }
    
    # 5. Store the filtered data in the list
    result_list[[result_counter]] <- filtered_data
    result_counter <- result_counter + 1
  }
}

# Combine all data frames in the list into a single data frame
final_result <- bind_rows(result_list)

# Print and save the result
print(final_result)
write.csv(final_result, "filtered_datasets.csv", row.names = FALSE)

#Verification
dataset_summary_exclusive <- final_result %>%
  group_by(dataset.id, Best.Annotation.Confidence.Level) %>%
  summarize(unique_smiles_count = n(), .groups = "drop")
print(dataset_summary_exclusive)
write.csv(dataset_summary_exclusive, "unique_smiles.csv", row.names = FALSE)

#Count USIs per dataset.id and confidence.level (no change here)
usi_counts <- datasets %>%
  group_by(dataset.id, Best.Annotation.Confidence.Level) %>%
  summarize(total_usi = n_distinct(feature.usi), .groups = "drop")

write_csv(usi_counts, "usi_counts.csv")
print(usi_counts)

#------------------------------------------------------------------------------#

##Compound Class Plot -- By Siyao Liu
# Count features per Best.Annotation.Compound.Class and dataset.id

# Read the filtered dataset and clean invalid superclasses
filtered_datasets <- final_result %>%
  filter(!(Best.Annotation.Compound.Class %in% c("None", "Others", "N/A", "NA", "")) & !is.na(Best.Annotation.Compound.Class))

plot_data <- filtered_datasets %>%
  group_by(Best.Annotation.Compound.Class, dataset.id) %>%
  summarise(Feature_Count = n(), .groups = "drop")

# Create and save the plot in one go
ggplot(plot_data, aes(x = Feature_Count, y = Best.Annotation.Compound.Class, fill = dataset.id)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "Unique Feature Count by Superclass and Dataset",
    x = "Number of Annotations",
    y = "Superclass",
    fill = "Dataset ID"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    legend.position = "right"
  ) -> p  # assign plot to object 'p'

# Save the plot
ggsave(
  filename = "Best_Annotation_Compound_Class_Annotations_counts_unique.svg",
  plot = p,
  width = 10,
  height = 14,
  dpi = 300
)


##Compound Class Plot -- By Siyao Liu
# Count features per Best.Annotation.Compound.Class and dataset.id

# Read the filtered dataset and clean invalid superclasses
datasets <- datasets %>%
  filter(!(Best.Annotation.Compound.Class %in% c("None", "Others", "N/A", "NA", "")) & !is.na(Best.Annotation.Compound.Class))

plot_data <- datasets %>%
  group_by(Best.Annotation.Compound.Class, dataset.id) %>%
  summarise(Feature_Count = n(), .groups = "drop")

# Create and save the plot in one go
ggplot(plot_data, aes(x = Feature_Count, y = Best.Annotation.Compound.Class, fill = dataset.id)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "Total Annotations by Superclass and Dataset",
    x = "Number of Annotations",
    y = "Superclass",
    fill = "Dataset ID"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    legend.position = "right"
  ) -> p  # assign plot to object 'p'

# Save the plot
ggsave(
  filename = "Best_Annotation_Compound_Class_Annotation_counts_total.svg",
  plot = p,
  width = 10,
  height = 14,
  dpi = 300
)