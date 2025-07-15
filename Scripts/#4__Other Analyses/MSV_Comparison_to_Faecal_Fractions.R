# Load necessary libraries
library(dplyr)
library(VennDiagram)
library(ggplot2) # Often useful for plotting

# Load your data
d128 <- read.csv("d128_post_standardisation.csv")
total <- read.csv("Y:/MA_BPA_Microbiome/Total-List-Of-Annotations.csv")

# --- MODIFICATION START ---

# 1. Handle NAs in Best.Annotation column for both dataframes
# Remove rows where Best.Annotation is NA
d128_cleaned <- d128 %>%
  filter(!is.na(Best.Annotation))

total_cleaned <- total %>%
  filter(!is.na(Best.Annotation))

# --- MODIFICATION END ---

# 2. Extract unique Best.Annotation values from the CLEANED dataframes
d128_annotations <- d128_cleaned %>%
  select(Best.Annotation) %>%
  distinct(Best.Annotation) %>%
  pull(Best.Annotation) # Convert to a vector

total_annotations <- total_cleaned %>%
  select(Best.Annotation) %>%
  distinct(Best.Annotation) %>%
  pull(Best.Annotation) # Convert to a vector

# 3. Prepare the list for VennDiagram
list_of_annotations <- list(
  d128 = d128_annotations,
  Total_List = total_annotations
)

# 4. Create the Venn Diagram
# Using the same parameters as in your provided code snippet
png("venn_diagram_annotations.png") # Start PNG device
venn.diagram(
  x = list_of_annotations,
  category.names = c("d128 Annotations", "Total Annotations"),
  filename = NULL, # Set to NULL because we're using the png() function directly
  output = TRUE,
  col = c("#4C0000", "#1A1A1A"), # Colors of the circles
  fill = c("#CC6666", "#333333"), # Fill colors of the circles
  alpha = 0.50, # Transparency of the fill
  cex = 1.5, # Size of set labels
  fontfamily = "sans",
  cat.cex = 1.5, # Size of category labels
  cat.fontfamily = "sans"
)
dev.off() # Close the PNG device

# Get the sizes for textual output (optional but good for verification)
num_d128 <- length(d128_annotations)
num_total <- length(total_annotations)
intersection <- length(intersect(d128_annotations, total_annotations))
only_d128 <- num_d128 - intersection
only_total <- num_total - intersection

cat("Number of unique annotations in d128 (NA-removed):", num_d128, "\n")
cat("Number of unique annotations in Total List (NA-removed):", num_total, "\n")
cat("Number of common annotations:", intersection, "\n")
cat("Annotations only in d128:", only_d128, "\n")
cat("Annotations only in Total List:", only_total, "\n")




### Part 2



d128 <- read.csv("Y:/MA_BPA_Microbiome/Dataset-Annotations/HGMD_0128.csv")

# Load necessary libraries
library(dplyr)

# Ensure Best.Annotation.Confidence.Level exists in your d128 dataframe.
# If it doesn't, you might have made a typo or need to create it.
# For example, if Best.Annotation.Confidence.Level is derived from Best.Annotation,
# you would need to define that logic first.
# Assuming 'Best.Annotation.Confidence.Level' is a direct column in d128.csv

# Calculate the number of entries for each Best.Annotation.Confidence.Level
annotation_level_summary <- d128 %>%
  # Handle potential NAs in Best.Annotation.Confidence.Level if you don't want them in the summary
  # filter(!is.na(Best.Annotation.Confidence.Level)) %>% # Uncomment this line if you want to exclude NAs
  
  group_by(Best.Annotation.Confidence.Level) %>%
  summarise(Count = n()) %>%
  ungroup() # Ungroup to remove the grouping attribute

# Print the summary
print(annotation_level_summary)

# You can also order it by count for better readability
annotation_level_summary_ordered <- annotation_level_summary %>%
  arrange(desc(Count))

print(annotation_level_summary_ordered)



###


# Filter for Best.Annotation.Confidence.Level == 2 and then count distinct Best.Annotation
level_2_distinct_annotations_count <- d128 %>%
  filter(Best.Annotation.Confidence.Level == 2) %>%
  # Optionally, filter out NAs in Best.Annotation if they exist and you don't want them counted
  # filter(!is.na(Best.Annotation)) %>%
  distinct(Best.Annotation) %>%
  nrow() # Count the number of rows (which are now distinct annotations)

# Print the result
cat("Number of distinct Best.Annotation entries where Best.Annotation.Confidence.Level is 2:",
    level_2_distinct_annotations_count, "\n")

# If you also want to see the distinct annotations themselves, you can do this:
level_2_distinct_annotations <- d128 %>%
  filter(Best.Annotation.Confidence.Level == 2) %>%
  # filter(!is.na(Best.Annotation)) %>% # Uncomment if needed
  distinct(Best.Annotation)

write.csv(level_2_distinct_annotations, "d128.csv")






## Part 3


# Remove entries where Best.Annotation contains "Candidate" (case-insensitive)
d128_filtered <- d128.2 %>%
  #filter(!str_detect(Best.Annotation, regex("Candidate", ignore_case = TRUE)))
  # Alternatively, using base R's grepl:
  filter(!grepl("Candidate", Best.Annotation, ignore.case = TRUE))

# You can check the dimensions before and after to see how many rows were removed
cat("Original number of rows in d128:", nrow(d128.2), "\n")
cat("Number of rows in d128 after removing 'Candidate' entries:", nrow(d128_filtered), "\n")
