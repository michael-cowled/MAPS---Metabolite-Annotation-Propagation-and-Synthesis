### META_COMPARE.R ###
##To compare metadata

library(tidyverse)
library(ggvenn)
library(dplyr)

# Specify the dataset of interest
dataset.id <- "HGMD_0125"

#metadata - GNPS compatible file (tab delimited .txt file)
metadata <- read_tsv("metadata.txt")

###---User-defined metadata categories---###

group.1 <- filter(metadata, ATTRIBUTE_Group == "Large-scale")
group.2 <- filter(metadata, ATTRIBUTE_Group == "Small-scale")
name.group.1 <- "Large-scale"
name.group.2 <- "Small-scale"

###--------------------------End of user input---------------------------###

final.annotation.df <- read.csv(paste0("Y:/MA_BPA_Microbiome/Dataset-Annotations/", dataset.id, ".csv"))

final.annotation.df3 <- final.annotation.df
final.annotation.df4 <- filter(final.annotation.df, !is.na(Best.Annotation))

# Initialize the new column with FALSE
final.annotation.df4$Containsgroup.1Filename <- FALSE

# Loop through each filename in group.1$filename
for (i in 1:nrow(group.1)) {
  filename <- group.1$filename[i] 
  print(filename)
  final.annotation.df4$Containsgroup.1Filename <- final.annotation.df4$Containsgroup.1Filename | grepl(filename, final.annotation.df4$Samples)
}

# Initialize the new column with FALSE
final.annotation.df4$Containsgroup.2Filename <- FALSE

# Loop through each filename in group.1$filename
for (i in 1:nrow(group.2)) {
  filename <- group.2$filename[i] 
  print(filename)
  final.annotation.df4$Containsgroup.2Filename <- final.annotation.df4$Containsgroup.2Filename | grepl(filename, final.annotation.df4$Samples)
}

#Create a column for presence in ONE (1) or BOTH (2):
final.annotation.df4 <- mutate(final.annotation.df4, Unique.1.OR.Mutual.2 = Containsgroup.1Filename + Containsgroup.2Filename)

Name <- c(name.group.1, name.group.2)

final.annotation.df4 %>%
  ggplot() +
  geom_venn(aes(A = Containsgroup.1Filename, B = Containsgroup.2Filename),
            set_names = Name,
            fill_color = c("#4C72B0", "#59A14F"), # More aesthetically pleasing colors
            fill_alpha = 0.7,                   # Add some transparency
            stroke_color = "grey30",            # Darker stroke for definition
            stroke_size = 0.8,                  # Slightly thicker stroke
            text_color = "white",
            text_size = 6,                      # Increased size for intersection numbers
            set_name_color = "black",           # Color for set names
            set_name_size = 7,                  # Increased size for set names
            show_percentage = TRUE,             # Show percentages within the diagram
            digits = 1) +                       # Number of decimal places for percentages
  theme_void() +
  labs(title = "Overlap of Group 1 and Group 2 Filenames") + # Add a title
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), # Center and style title
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm") # Add a bit more margin around the plot
  )



df_group1 <- final.annotation.df4 %>%
  filter(Containsgroup.1Filename == TRUE)

# Data frame for annotations related to Small-scale samples
df_group2 <- final.annotation.df4 %>%
  filter(Containsgroup.2Filename == TRUE)

# Define the desired column order
desired_cols <- c(
  "dataset.id",
  "feature.usi",
  "Samples",
  "Best.Annotation",
  "Best.Annotation.Smiles",
  "Best.Annotation.Confidence.Level",
  "Best.Annotation.Compound.Class"
)

df_group1_tidy <- df_group1 %>%
  mutate(dataset.id = name.group.1) %>%
  select(all_of(desired_cols))

df_group2_tidy <- df_group2%>%
  mutate(dataset.id = name.group.2) %>%
  select(all_of(desired_cols))

# 3. Combine the two dataframes
datasets <- bind_rows(
  df_group1_tidy,
  df_group2_tidy
)

datasets$Best.Annotation.Confidence.Level <- as.numeric(as.character(datasets$Best.Annotation.Confidence.Level))

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
## Compound Class Plot -- By Siyao Liu
# Count features per Best.Annotation.Compound.Class and dataset.id

# Read the filtered dataset and clean invalid superclasses
filtered_datasets <- final_result %>%
  filter(!(Best.Annotation.Compound.Class %in% c("None", "Others", "N/A", "NA", "")) & !is.na(Best.Annotation.Compound.Class))

plot_data <- filtered_datasets %>%
  group_by(Best.Annotation.Compound.Class, dataset.id) %>%
  summarise(Feature_Count = n(), .groups = "drop")

# Create and save the plot in one go     ###UNIQUE
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

# Create and save the plot in one go      ###TOTAL
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

print(p)