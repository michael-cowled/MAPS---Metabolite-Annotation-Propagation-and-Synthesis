#Making Proportion of superclasses per apple figure

library(ggplot2)
library(reshape2) 
library(dplyr)
library(tidyr)

#Read in CSV's and merge
class.data <- read.csv("ms2query.csv")                                               #####!!!###### Check spelling of filename
names(class.data)[1] <- "row.ID"
class.data <- select(class.data, "row.ID", "npc_superclass_results")

quant.data <- read.csv("quant.csv")                                      #####!!!######  Check spelling of filename
quant.data <- quant.data[, -c(2:13)] # Removes unneccessary columns

#Tidy filenames:
# Get the column names of quant.data
colnames_quant <- colnames(quant.data)

# Remove the suffix ".mzML.Peak.area" if present
colnames_quant <- sub("\\.mzML\\.Peak\\.area$", "", colnames_quant)

# Assign the modified column names back to quant.data
colnames(quant.data) <- colnames_quant

#Merge data
merged_data <- merge(class.data, quant.data, by = "row.ID") %>%
  select(-row.ID)

#Tidy data to be presence/absence based for feature. Edit threshold value to be considered a feature?
# Select columns 3 to last column
cols_to_replace <- 2:ncol(merged_data)

# Replace values greater than 0 with 1
merged_data[, cols_to_replace] <- apply(merged_data[, cols_to_replace], 2, function(x) ifelse(x > 0, 1, x))


#Function: SUMMARISE_SUPERCLASSES:
summarise_superclasses <- function(j, merged_data) {
  superclasses <- unique(merged_data$npc_superclass_results)
  temp_data_frame <- data.frame()  # Initialize an empty dataframe to store results
  
  # Iterate through each unique value in superclasses
  for (i in superclasses) {
    print(i)
    # Filter merged_data based on npc_superclass_results == i
    filt_data <- merged_data[merged_data$npc_superclass_results == i, c(1, j)]
    
    # Calculate sum of the filtered data
    total_sum <- sum(filt_data[, 2], na.rm = TRUE)  # Assuming column 2 contains numeric values
    
    # Create a temporary dataframe to store the results for this iteration
    temp_df <- data.frame(npc_superclass_results = i, total_sum = total_sum)
    
    # Rename colnames
    colnames(temp_df) <- c("npc_superclass_results", colnames(merged_data)[j])
    
    # Append the temporary dataframe to the main dataframe
    temp_data_frame <- rbind(temp_data_frame, temp_df)
  }
  
  return(temp_data_frame)
}

#Now, to generate for all samples:
# Initialize the merged dataframe with the result from the first iteration
merged_df <- summarise_superclasses(j = 2, merged_data)

# Loop over columns 3 to ncol(merged_data)
for (j in 3:ncol(merged_data)) {
  # Call the summarise_superclasses function for each column
  tdf <- summarise_superclasses(j, merged_data)
  
  # Merge the current dataframe with the merged dataframe
  merged_df <- merge(merged_df, tdf, by = "npc_superclass_results", all = TRUE)
}

# Create an empty dataframe with the same structure as merged_df
merged_df_perc <- data.frame(matrix(NA, nrow = nrow(merged_df), ncol = ncol(merged_df)))
colnames(merged_df_perc) <- colnames(merged_df)

# Loop over each column in merged_df
for (col in names(merged_df)[-1]) {
  # Calculate the sum of values in the column
  col_sum <- sum(merged_df[[col]], na.rm = TRUE)
  
  # Replace each value in the column with the percentage
  merged_df_perc[[col]] <- (merged_df[[col]] / col_sum) * 100
}

merged_df_perc[,1] <- merged_df[,1]

## Change from wide to long df format
library(tidyr)

long <- merged_df_perc %>% 
  pivot_longer(
    cols = -npc_superclass_results,  # Exclude the 'npc_superclass_results' column
    names_to = "Sample",
    values_to = "Perc"
  )


#Generate bubble chart figure

category_colors <- rainbow(length(unique(long$'npc_superclass_results')))

# Create the plot object with rotated x-axis labels
plot_obj <- ggplot(long, aes(x = Sample, y = npc_superclass_results, size = Perc, color = npc_superclass_results)) +
  geom_point() +
  scale_size(range = c(1, 8)) +
  scale_color_manual(values = category_colors) +
  labs(x = "Filename", y = "Superclasses", size = "Proportion of Metabolites") +
  guides(color = FALSE) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x-axis labels
  scale_y_discrete(limits = rev(levels(long$npc_superclass_results)))

# Print the plot
print(plot_obj)





###OPTIONAL CODE: filters be a minimum of 1 to be shown and present

library(dplyr)

# Filter the data to remove rows where Perc is not more than 0
filtered_long <- long %>% 
  filter(Perc > 0)

# Create the plot object with modified x-axis labels
plot_obj <- ggplot(filtered_long, aes(x = Sample, y = npc_superclass_results, size = Perc, color = npc_superclass_results)) +
  geom_point() +
  scale_size(range = c(1, 8)) +
  scale_color_manual(values = category_colors) +
  labs(x = "Filename", y = "Superclasses", size = "Proportion of Metabolites") +
  guides(color = FALSE) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x-axis labels
  scale_y_discrete(limits = rev(levels(filtered_long$npc_superclass_results)))

# Print the plot
print(plot_obj)
