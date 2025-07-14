#BUBBLE CHART: FOR MS2QUERY                                
library(ggplot2)
library(reshape2) 
library(dplyr)
library(tidyr)
library(readr)
library(svglite)
library(readxl)

##---##

sheet_names <- excel_sheets(excel_file) # Get the sheet names
for (sheet in sheet_names) {
  data <- read_excel(excel_file, sheet = sheet)
  write.csv(data, file = paste0(sheet, ".csv"), row.names = FALSE)
}

metadata <- read.csv("A - Analysis.csv")
names(metadata)[1] <- "Sample"
metadata <- select(metadata, Sample, Book.Code)

long <- long %>%
  left_join(metadata, by = "Sample") %>%
  mutate(Sample = coalesce(Book.Code, Sample)) %>%
  select(-Book.Code)

#Read in CSV's and merge
class.data <- read.csv("ms2query.csv") %>%                            #####!!!###### Check spelling of filename
  filter(ms2query_model_prediction >= 0.63)
names(class.data)[10] <- "row.ID"
class.data <- select(class.data, "row.ID", "npc_superclass_results") %>%
  filter(npc_superclass_results != "None")

quant.data <- read.csv("DATA_iimn_gnps_quant.csv")                    #####!!!######  Check spelling of filename
quant.data <- quant.data[, -c(2:13)] # Removes unneccessary columns

quant.data[, 2:ncol(quant.data)][quant.data[, 2:ncol(quant.data)] <10000] <- 0   ###SELECT NOISE THRESHOLD: Note: this is AREA, not height

#Tidy filenames:
class.data <- class.data %>%
  mutate(npc_superclass_results = ifelse(npc_superclass_results == "", "Unclassified", npc_superclass_results)) %>%
  filter(npc_superclass_results != "Unclassified")

# Get the column names of quant.data
colnames_quant <- colnames(quant.data)

# Remove the suffix ".mzML.Peak.area" if present
colnames_quant <- sub("\\.mzML\\.Peak\\.area$", "", colnames_quant)

# Assign the modified column names back to quant.data
colnames(quant.data) <- colnames_quant

#Merge data 
QUERY_TABLE <- merge(class.data, quant.data, by = "row.ID") #%>%                 ##Captures a table to look at FEATURES
#  filter(Prerun_Blank_02 <= 10 | Midrun_Blank_01 <= 10, Postrun_Blank_01 <= 10)    ##Change BLANK NAMES here! -- Removed features in the blanks
#filter(Prerun_blank_01 <= 10)

merged_data <- select(QUERY_TABLE, -row.ID) #%>%   
#  select(-Prerun_Blank_02, -Midrun_Blank_01, -Postrun_Blank_01)                  ##Now remove blanks

##EXAMPLE QUERY: 
# query <- QUERY_TABLE %>%
#select(npc_superclass_results, row.ID, HGMA_0414) %>% 
#        filter(npc_superclass_results == "Macrolides")
# View(query)

#Function: SUMMARISE_superclasss:
summarise_superclasss <- function(j, merged_data) {
  superclasss <- unique(merged_data$npc_superclass_results)
  temp_data_frame <- data.frame()  # Initialize an empty dataframe to store results
  
  # Iterate through each unique value in superclasss
  for (i in superclasss) {
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
merged_df <- summarise_superclasss(j = 2, merged_data)

# Loop over columns 3 to ncol(merged_data)
for (j in 3:ncol(merged_data)) {
  # Call the summarise_superclasss function for each column
  tdf <- summarise_superclasss(j, merged_data)
  
  # Merge the current dataframe with the merged dataframe
  merged_df <- merge(merged_df, tdf, by = "npc_superclass_results", all = TRUE)
}

## Change from wide to long df format 
library(tidyr)

long <- merged_df %>% 
  pivot_longer(
    cols = -npc_superclass_results,  # Exclude the 'npc_superclass_results' column
    names_to = "Sample",
    values_to = "Perc"
  )

#Generate bubble chart figure

category_colors <- rainbow(length(unique(long$'npc_superclass_results')))

library(dplyr)

# Filter the data to remove rows where Perc is not more than 0
filtered_long <- long %>% 
  filter(Perc > 0)

# Create the plot object with modified x-axis labels
plot_obj <- ggplot(filtered_long, aes(x = Sample, y = npc_superclass_results, size = Perc, color = npc_superclass_results)) +
  geom_point() +
  scale_size(range = c(1, 8)) +
  scale_color_manual(values = category_colors) +
  labs(x = "Fraction", y = "Metabolite Class", size = "Sum of Peak Areas") +
  guides(color = "none") +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),  # Increase x-axis text size
    axis.text.y = element_text(size = 12),  # Increase y-axis text size
    axis.title.x = element_text(size = 14),  # Increase x-axis title size
    axis.title.y = element_text(size = 14),  # Increase y-axis title size
    legend.text = element_text(size = 12),  # Increase legend text size
    plot.title = element_text(size = 16)  # Increase plot title size if used
  ) +
  scale_y_discrete(limits = rev(levels(filtered_long$npc_superclass_results)))

# Print the plot
print(plot_obj)

ggsave(filename = "bubblechart.svg", 
       plot = plot_obj, 
       width = 50, 
       height = 40, 
       units = "cm")