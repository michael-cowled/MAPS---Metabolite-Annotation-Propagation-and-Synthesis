### Cumulative-histogram.R
##To make a cumulative histogram to show unique annotations with respect to crude and previous fractions

# Define the crude sample (or samples) to act as the baseline                     ###USER-INPUT###
remove_sample <- "HGMA_1598"

# Filter out rows containing the specified sample
Annotations.with.samples <- Annotations.with.samples %>%
  filter(!str_detect(Samples, remove_sample))

library(tidyverse)

# Initialize an empty dataframe to store the results
cumulative_annotations <- data.frame()

# Get a sorted list of unique sample numbers
sample_numbers <- sort(unique(as.numeric(str_extract(unlist(str_split(Annotations.with.samples$Samples, "; ")), "\\d+"))))

# Loop through each sample number
for (i in seq_along(sample_numbers)) {
  # Filter data for the current sample number
  current_sample_data <- Annotations.with.samples %>%
    filter(str_detect(Samples, paste0("_", sample_numbers[i])))
  
  # Calculate the number of new annotations in the current sample
  new_annotations <- sum(current_sample_data$annotation)
  
  # Add the new annotations to the cumulative count, even if 0
  if (nrow(cumulative_annotations) == 0) {
    cumulative_annotations <- data.frame(
      Sample.Number = i,
      Cumulative.Unique.Annotations = new_annotations
    )
  } else {
    cumulative_annotations <- cumulative_annotations %>%
      add_row(
        Sample.Number = i,
        Cumulative.Unique.Annotations = new_annotations + last(cumulative_annotations$Cumulative.Unique.Annotations)
      )
  }
  
  # Remove the processed features from the main dataframe
  Annotations.with.samples <- Annotations.with.samples %>%
    filter(!feature.ID %in% current_sample_data$feature.ID)
}

# Ensure all chronological sample numbers are present
all_sample_numbers <- data.frame(Sample.Number = seq_along(sample_numbers))
cumulative_annotations <- full_join(all_sample_numbers, cumulative_annotations, by = "Sample.Number")

# Fill missing Cumulative.Unique.Annotations with the previous value
cumulative_annotations <- cumulative_annotations %>%
  fill(Cumulative.Unique.Annotations, .direction = "down")

# Print the result
print(cumulative_annotations)

# Create the cumulative histogram
plot_obj3 <- ggplot(cumulative_annotations, aes(x = Sample.Number, y = Cumulative.Unique.Annotations)) +
  geom_col() +
  labs(
    x = "Fraction Number",
    y = "Cumulative Number of Unique Annotations"
  ) +
  theme(text = element_text(size = 24))

print(plot_obj3)

ggsave(filename = "cumulative-histogram.svg", 
       plot = plot_obj3, 
       width = 50, 
       height = 40, 
       units = "cm")
