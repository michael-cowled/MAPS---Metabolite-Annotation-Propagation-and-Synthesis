library(ggplot2)
library(dplyr)
library(readr)
#---------USER GUIDE: SET THE DIRECTORY TO THE DATASET COMPARIOSON FOLDER, AND CHANGE THE DATASET NAME ACCORDINGLY IN LINE 8 - 13, AND LINE 17-18----------#
# Set path
base_path <- "Y:/MA_BPA_Microbiome/Dataset-Annotations"

# Load and label datasets
df_0110 <- read_csv(file.path(base_path, "HGMD_0112.csv")) %>%
  mutate(Dataset = "HGMD_0112")

df_0111 <- read_csv(file.path(base_path, "HGMD_0113.csv")) %>%
  mutate(Dataset = "HGMD_0113")



# Combine datasets
combined_df <- bind_rows(df_0110, df_0111)


# Define a low-saturation (pastel) color palette
custom_colors_20 <- c(
  "#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854",
  "#ffd92f", "#e5c494", "#b3b3b3", "#a6cee3", "#1f78b4",
  "#b2df8a", "#33a02c", "#fb9a99", "#bb4b3d", "#fdbf6f",
  "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928", "#A6CEE3", "#B2DF8A", "#FDBF6F", "#CAB2D6", "#FB9A99"
)

##------------------------Stacked density plot---------------------------##
# Plot: vertically stacked facets with shared x-axis
# Find the maximum rt
max_rt <- max(combined_df$rt, na.rm = TRUE)

# Plot with manual x-axis limit
p1 <- ggplot(combined_df, aes(x = rt, fill = Dataset)) +
  geom_density(alpha = 0.3, adjust = 1.2, color = NA) +
  facet_wrap(~Dataset, ncol = 1, scales = "fixed") +
  labs(title = "Feature Density along Retention Time",
       x = "Retention Time (rt)",
       y = "Density") +
  scale_fill_brewer(palette = "Set1") +
  coord_cartesian(xlim = c(0, max_rt)) +  # restrict x axis to data range
  theme_minimal() +
  theme(strip.text = element_text(size = 12, face = "bold"),
        legend.position = "none")
print(p1)

#-----------------------stacked bar chart showing superclass composition------------#
# Step 1: Remove unwanted classes
cleaned_df <- combined_df %>%
  filter(!(Best.Annotation.Compound.Class %in% c("None", "Others", "N/A", "NA", "")) & !is.na(Best.Annotation.Compound.Class))

# Step 2: Find top 15 classes per dataset after filtering
top_classes <- cleaned_df %>%
  group_by(Dataset, Best.Annotation.Compound.Class) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Dataset) %>%
  slice_max(order_by = Count, n = 15) %>%
  ungroup()

# Step 3: Filter cleaned_df to include only top classes
filtered_df <- inner_join(cleaned_df, top_classes, 
                          by = c("Dataset", "Best.Annotation.Compound.Class"))

# Compute counts, percentages, and angles
class_counts <- filtered_df %>%
  group_by(Dataset, Best.Annotation.Compound.Class) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Dataset) %>%
  mutate(Percent = Count / sum(Count) * 100,
         Label = if_else(row_number(desc(Count)) <= 5,
                         paste0(Best.Annotation.Compound.Class, "\n", round(Percent, 1), "%"),
                         NA_character_),
         ymax = cumsum(Count),
         ymin = lag(ymax, default = 0),
         mid = (ymax + ymin) / 2)


# Step 1: Reorder compound class factor by total count
class_order <- class_counts %>%
  group_by(Best.Annotation.Compound.Class) %>%
  summarise(Total = sum(Count), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  pull(Best.Annotation.Compound.Class)

class_counts$Best.Annotation.Compound.Class <- factor(
  class_counts$Best.Annotation.Compound.Class,
  levels = class_order
)

# Step 2: Plot with sorted legend
p2<-ggplot(class_counts, aes(x = Dataset, y = Count, fill = Best.Annotation.Compound.Class)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = custom_colors_20) +
  labs(title = "Feature Count by compound Class",
       x = "Dataset", y = "Number of Features",
       fill = "Compound Class") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

#----------------------Annotation confidence conparison------------------------#
# Filter and clean
confidence_data <- combined_df %>%
  filter(!Best.Annotation.Confidence.Level %in% c("NA", "N/A", "", "Unknown")) %>%
  mutate(Best.Annotation.Confidence.Level = factor(Best.Annotation.Confidence.Level, ordered = TRUE))

# Plot: Bar chart of confidence level counts per dataset
p3<-ggplot(confidence_data, aes(x = Best.Annotation.Confidence.Level, fill = Dataset)) +
  geom_bar(position = position_dodge(width = 0.9)) +
  geom_text(stat = "count",
            aes(label = ..count..),
            position = position_dodge(width = 0.9),
            vjust = -0.3,
            size = 3) +
  scale_fill_manual(values = custom_colors_20) +
  labs(title = "Annotation Confidence Level by Dataset",
       x = "Confidence Level", y = "Feature Count",
       fill = "Dataset") +
  theme_minimal(base_size = 12)

# Filter and prepare data
df_filtered <- combined_df %>%
  filter(!is.na(Best.Annotation.Compound.Class),
         !Best.Annotation.Compound.Class %in% c("", "None", "NA", "N/A"))

###-----------------------stacked bar plot showing the annotation source----##

percent_data <- df_filtered %>%
  count(Best.Annotation.Compound.Class, Best.Annotation.Type) %>%
  group_by(Best.Annotation.Compound.Class) %>%
  mutate(Percentage = n / sum(n) * 100) %>%
  ungroup()

stacked_bar_plot <- ggplot(percent_data, 
                           aes(x = reorder(Best.Annotation.Compound.Class, -Percentage), 
                               y = Percentage, 
                               fill = Best.Annotation.Type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values = custom_colors_20) +
  labs(title = "Percentage Contribution of Annotation Sources",
       subtitle = "By compound class",
       x = "Compound Class",
       y = "Percentage (%)",
       fill = "Annotation Source") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size= 6.5),
        plot.title = element_text(face = "bold", size = 14),
        legend.position = "right") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

# Display plots

print(stacked_bar_plot)

# Save plots (optional)
ggsave("stacked_bar_percentage.png", stacked_bar_plot, width = 10, height = 6, dpi = 600)

###------Ridge plot-----------####

# Bin RT into exact 1-minute intervals and count features
rt_binned <- df_filtered %>%
  mutate(
    RT_bin = floor(rt),  # Creates integer bins (0,1,2,...)
    RT_bin_label = sprintf("%d-%.4f", RT_bin, RT_bin + 0.9999)  # Precise labels
  ) %>%
  count(RT_bin, RT_bin_label, Best.Annotation.Compound.Class, name = "Count") %>%
  mutate(RT_mid = RT_bin + 0.5) %>%  # Midpoint for plotting
  group_by(Best.Annotation.Compound.Class) %>%
  mutate(Total_Class_Count = sum(Count)) %>%  # For ordering
  ungroup()

# Bubble plot visualization
rt.count <- max(rt_binned$RT_bin) + 1
bubble<-ggplot(rt_binned, 
               aes(x = RT_mid, 
                   y = reorder(Best.Annotation.Compound.Class, Total_Class_Count),
                   size = Count, 
                   color = Best.Annotation.Compound.Class)) +
  geom_point(alpha = 0.8) +
  scale_x_continuous(breaks = 0:rt.count, limits = c(0, rt.count)) +
  scale_size_continuous(range = c(3, 15), name = "Feature Count") +
  scale_color_discrete() +
  labs(title = "Retention Time Distribution by Compound Class",
       x = "Retention Time (minutes)",
       y = "Compound Class") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5),
        panel.grid.major.y = element_line(color = "grey90")) +
  guides(color = "none") # This line removes the legend for the 'color' aesthetic
print(bubble)
# Save plots (optional)
ggsave("rt distribution of compound class.png", bubble, width = 20, height = 10, dpi = 600)

#-------------------------save the plot--Depending on what plot is needed-------------------------------#
library(patchwork)
library(ggplot2)
combined_plot <-  (p1)/(p2 | p3)
ggsave("combined_figure.png", combined_plot, width = 14, height = 12, dpi = 600)