#BUBBLE CHART: FOR SIRIUS                                
library(ggplot2)
library(reshape2) 
library(dplyr)
library(tidyr)
library(readr)

#Read in CSV's and merge
class.data <-readr::read_tsv("canopus_structure_summary.tsv")                      #####!!!###### Check spelling of filename
names(class.data)[c(7, 8,26)] <- c("npc_superclass_results", "prob", "row.ID")
class.data <- filter(class.data, prob >= 0.6)
class.data <- select(class.data, "row.ID", "npc_superclass_results")

quant.data <- read.csv("DATA_iimn_gnps_quant.csv")                    #####!!!######  Check spelling of filename
quant.data <- quant.data[, -c(2:13)] # Removes unneccessary column

quant.data[, 2:ncol(quant.data)][quant.data[, 2:ncol(quant.data)] < 1000] <- 0   ###SELECT NOISE THRESHOLD

#Tidy filenames:
class.data <- class.data %>%
  mutate(npc_superclass_results = ifelse(npc_superclass_results == "", "Unclassified", npc_superclass_results))

# Get the column names of quant.data
colnames_quant <- colnames(quant.data)

# Remove the suffix ".mzML.Peak.area" if present
colnames_quant <- sub("\\.mzML\\.Peak\\.area$", "", colnames_quant)

# Assign the modified column names back to quant.data
colnames(quant.data) <- colnames_quant

#Merge data 
QUERY_TABLE <- merge(class.data, quant.data, by = "row.ID") %>%                 ##Captures a table to look at FEATURES
  filter(Prerun_blank_01 <= 10 | Prerun_blank_02 <= 10 | Prerun_blank_03 <= 10 | Prerun_blank_04 <= 10)    ##Change BLANK NAMES here! -- Removed features in the blanks
  #filter(Prerun_blank_02 <= 10)

merged_data <- select(QUERY_TABLE, -row.ID) %>%   
  select(-Prerun_blank_01, -Prerun_blank_02, -Prerun_blank_03, -Prerun_blank_04)                  ##NOw remove blanks

##EXAMPLE QUERY: 
#query <- select(QUERY_TABLE, HGMA_0381, npc_superclass_results, row.ID) %>% 
#          filter(npc_superclass_results == "Fatty Acids and Conjugates", HGMA_0381 > 0)
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


###Mapping to metadata                                                            ###Download latest version as CSV,
###And ensure codes have -01, -02 format
metadata <- read.csv("HGM - Data Management System.csv")
names(metadata)[1] <- "Sample"
metadata <- select(metadata, Sample, Book.Code)

long <- long %>%
  left_join(metadata, by = "Sample") %>%
  mutate(Sample = coalesce(Book.Code, Sample)) %>%
  select(-Book.Code)

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
  labs(x = "Fraction", y = "Superclass", size = "Cum. Sum of Peak Areas") +
  guides(color = "none") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +  # Rotate x-axis labels
  scale_y_discrete(limits = rev(levels(filtered_long$npc_superclass_results)))

# Print the plot
print(plot_obj)
