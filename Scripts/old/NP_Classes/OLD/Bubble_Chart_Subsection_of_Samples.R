#AFTER 'long' is generated

# Define the pattern to match "-15" to "-17"
pattern <- "-1[5-7]"

# Filter the dataframe
long_filtered <- long %>%
  filter(grepl(pattern, Sample))

# Print the filtered dataframe
print(long_filtered)

#Generate bubble chart figure

category_colors <- rainbow(length(unique(long_filtered$'npc_superclass_results')))

library(dplyr)

# Filter the data to remove rows where Perc is not more than 0
filtered_long <- long_filtered %>% 
  filter(Perc > 0)

# Create the plot object with modified x-axis labels
plot_obj <- ggplot(filtered_long, aes(x = Sample, y = npc_superclass_results, size = Perc, color = npc_superclass_results)) +
  geom_point() +
  scale_size(range = c(1, 8)) +
  scale_color_manual(values = category_colors) +
  labs(x = "Filename", y = "superclasss", size = "Cum. Sum of Peak Areas") +
  guides(color = "none") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +  # Rotate x-axis labels
  scale_y_discrete(limits = rev(levels(filtered_long$npc_superclass_results)))

# Print the plot
print(plot_obj)
