# Install and load necessary package if you haven't already
# install.packages("ggplot2")
# install.packages("dplyr") # For data manipulation

library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(ComplexUpset)
library(svglite)
library(tidyverse)

#Using df 'combined' resulting from Extract annotations based on metadata.R
##Remember to standardise and deduplicate beforehand (by method)

annotation_cleaned <- combined %>%
  filter(confidence.level %in% 1:5) %>%
  mutate(
    Method = case_when(
      str_detect(metadata, "HILIC") ~ "HILIC",
      str_detect(metadata, "Phe-Hex") ~ "Phe-Hex"
    ),
    Scale = case_when(
      str_detect(metadata, "_LS") ~ "Large-Scale",
      str_detect(metadata, "_SS") ~ "Small-Scale"
    ),
    Confidence_Level = confidence.level
  ) %>%
  count(Confidence_Level, Method, Scale, name = "Number_of_Compounds")

# Step 2: Ensure all combinations exist, filling with 0 if missing
data_df <- annotation_cleaned %>%
  complete(
    Confidence_Level = 1:5,
    Method = c("Phe-Hex", "HILIC"),
    Scale = c("Small-Scale", "Large-Scale"),
    fill = list(Number_of_Compounds = 0)
  ) %>%
  arrange(Confidence_Level, Method, Scale)

# Filter out Confidence_Level 1
data_df <- data_df %>%
  filter(Confidence_Level != 1)

# Convert Confidence_Level to a factor for proper x-axis labeling
data_df$Confidence_Level <- as.factor(data_df$Confidence_Level) %>%
  droplevels()


# Create a new grouping variable for the fill aesthetic
data_df <- data_df %>%
  mutate(Bar_Group = paste(Method, Scale, sep = "_")) %>%
  # Ensure the desired order for bars within each Confidence Level group
  mutate(Bar_Group = factor(Bar_Group, levels = c(
    "Phe-Hex_Small-Scale", "Phe-Hex_Large-Scale",
    "HILIC_Small-Scale", "HILIC_Large-Scale"
  )))


# Define custom colors
colors_manual <- c(
  "Phe-Hex_Small-Scale" = "#ADD8E6", # Light Blue
  "Phe-Hex_Large-Scale" = "#4682B4", # Dark Blue (Steel Blue)
  "HILIC_Small-Scale" = "#90EE90",   # Light Green
  "HILIC_Large-Scale" = "#228B22"    # Dark Green (Forest Green)
)

# 2. Create the grouped bar chart
# Create the grouped bar chart without labels on top of the bars
# Create the grouped bar chart with a y-axis limit
ggplot(data_df, aes(x = Confidence_Level, y = Number_of_Compounds, fill = Bar_Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_manual(
    values = colors_manual,
    name = "Method and Scale", # Legend title
    labels = c(
      "Phe-Hex_Small-Scale" = "Phe-Hex (Small-Scale)",
      "Phe-Hex_Large-Scale" = "Phe-Hex (Large-Scale)",
      "HILIC_Small-Scale" = "HILIC (Small-Scale)",
      "HILIC_Large-Scale" = "HILIC (Large-Scale)"
    )
  ) +
  labs(
    title = "Comparison of Compound Counts by Confidence Level, Method, and Scale (Excluding Level 1)",
    x = "Confidence Level",
    y = "Number of Compounds"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 32),
    axis.text.x = element_text(size = 24),
    axis.text.y = element_text(size = 24),
    axis.title = element_text(size = 28, face = "bold"),
    legend.title = element_text(size = 26, face = "bold"),
    legend.text = element_text(size = 24),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  coord_cartesian(ylim = c(0, 10000))

# IMPORTANT: You will almost certainly need to increase the width and height
# when saving the plot, otherwise text will overlap significantly.
# A good starting point would be to double the width/height as well.
ggsave("compound_counts_grouped_bar_chart_doubled_fonts.svg", width = 16, height = 10, dpi = 300)

# Compound Class Plot (Total Features) - from 'all_combined_datasets' before filtering for uniqueness
plot_data_total <- combined %>%
  filter(!(NPC.superclass %in% c("None", "Others", "N/A", "NA", "")) & !is.na(NPC.superclass)) %>%
  group_by(NPC.superclass, metadata) %>%
  summarise(Feature_Count = n(), .groups = "drop") %>%
  # Filter categories with less than 10 annotations
  filter(Feature_Count >= 10)

ggplot(plot_data_total, aes(x = Feature_Count, y = NPC.superclass, fill = metadata)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "Total Annotations by Superclass and Dataset (All Four, >= 3 Annotations)",
    x = "Number of Annotations",
    y = "Superclass",
    fill = "metadata"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    legend.position = "right"
  ) -> p_total

ggsave(
  filename = "Compound_Class_Annotation_counts_total_all_four.svg",
  plot = p_total,
  width = 10,
  height = 14,
  dpi = 300
)


##Creation of Upset Plot
# Filter to non-NA SMILES
smiles_by_group <- combined %>%
  filter(!is.na(smiles)) %>%
  group_by(metadata) %>%
  summarise(smiles = list(unique(smiles))) %>%
  deframe()

# Get metadata groups (up to 4 for Venn)
groups <- names(smiles_by_group)

# Example: Compute pairwise overlaps and percent overlap matrix
overlap_matrix <- matrix(NA, nrow = length(groups), ncol = length(groups),
                         dimnames = list(groups, groups))

for (i in groups) {
  for (j in groups) {
    intersect_count <- length(intersect(smiles_by_group[[i]], smiles_by_group[[j]]))
    union_count <- length(union(smiles_by_group[[i]], smiles_by_group[[j]]))
    overlap_matrix[i, j] <- round(intersect_count / union_count * 100, 1)
  }
}

library(ComplexUpset)

# Make long format
upset_data <- combined %>%
  filter(!is.na(smiles)) %>%
  distinct(smiles, metadata) %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = metadata, values_from = present, values_fill = 0)

# Plot UpSet
upset(upset_data, intersect = names(upset_data)[-1])