####Plots based on Annotation-table####

##Annotation-table.R must be run first before generating plot.

              ###Feature-Counter.R counts the frequency of annotations (including propagated Ms/MS) across fractions and creates a barchart
              
              # Load necessary libraries
              library(dplyr)
              library(tidyr)
              library(ggplot2)
              
              ### Feature and Annotation counter per Sample/Fraction
              ## This follows up on Annotation-table.R and Cytoscape-tidier.R
              
              # 1. Preparing Annotations.with.samples
              Samples <- select(final.annotation.df, feature.ID, Samples, compound.name)  # Extract relevant sample data
              names(Samples)[names(Samples) == "compound.name"] <- "annotation"
              
              # Merge annotations with samples (before propagation)
              Annotations.with.samples <- Samples %>%
                mutate(annotation = ifelse(is.na(annotation), 0, 1))  # Convert 'annotation' column to presence/absence (0/1)
              
              # 2. Summarising annotations per sample (before propagation)
              
              # Sum of all annotations per sample
              sum_annotations_per_sample <- Annotations.with.samples %>%
                separate_rows(Samples, sep = "; ") %>%  # Split 'Samples' into individual sample rows
                group_by(Samples) %>%                   # Group by individual samples
                summarise(sum_annotations = sum(annotation, na.rm = TRUE))  # Sum annotations per sample
              
              # Sum of unique annotations per sample (where annotation is only for one sample)
              sum_unique_annotations_per_sample <- Annotations.with.samples %>%
                mutate(Samples_count = str_count(Samples, ";") + 1) %>%  # Count number of samples per annotation
                filter(Samples_count == 1) %>%                          # Filter to only unique annotations
                separate_rows(Samples, sep = "; ") %>%                  # Split samples into separate rows
                group_by(Samples) %>%                                   # Group by individual samples
                summarise(sum_unique_annotations = sum(annotation, na.rm = TRUE))  # Sum unique annotations per sample
              
              # 3. Merging all dataframes by 'Samples'
              
              # Merge dataframes on 'Samples'
              merged_df <- sum_annotations_per_sample %>%
                full_join(sum_unique_annotations_per_sample, by = "Samples")
              
              # 4. Adding total number of annotations per sample (counting both 1 and 0)
              
              # Calculate total number of features per sample (both 1 and 0) and add as a new column
              total_features_per_sample <- Annotations.with.samples %>%
                separate_rows(Samples, sep = "; ") %>%  # Split 'Samples' into individual sample rows
                group_by(Samples) %>%
                summarise(total_features = n())  # Count all annotations (both 1 and 0) per sample
              
              # Merge the total annotations column into the final dataframe
              merged_df <- merged_df %>%
                full_join(total_features_per_sample, by = "Samples")
              
              # View the final merged dataframe
              merged_df
              
              ###Mapping to metadata                                                            ###Download latest version as CSV,
              ###And ensure codes have -01, -02 format
              metadata <- read.csv("HGM/A - Analysis.csv")
              names(metadata)[1] <- "Sample"
              names(merged_df)[1] <- "Sample"
              metadata <- select(metadata, Sample, Book.Code)
              
              merged_df <- merged_df %>%
                full_join(metadata, by = "Sample") %>%
                mutate(Sample = coalesce(Book.Code, Sample)) %>%
                select(-Book.Code) %>%
                filter(!is.na(sum_annotations))
              
              write.csv(merged_df, paste0(folder, "/features-per-sample.csv"))
              
              #5. Barchart
              
              # Generate rainbow colors with alpha (transparency) set to 0.7
              category_colors <- rainbow(length(unique(merged_df$Sample)), alpha = 0.7)
              
              # Create the bar chart
              plot_obj2 <- ggplot(merged_df, aes(x = Sample, y = sum_annotations, fill = Sample)) +
                geom_col() +
                xlab("Sample") +
                ylab("Sum of Annotations") +
                scale_fill_manual(values = category_colors) +
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),  # Rotate x-axis labels
                      text = element_text(size = 24)) +  # Increase overall font size
                        ylim(0, 1000)  # Set y-axis limits
                        
              print(plot_obj2)

              ggsave(filename = paste0(folder, "/feature-barchart.svg"), 
                     plot = plot_obj2, 
                     width = 50, 
                     height = 40, 
                     units = "cm")

### Cumulative-histogram.R
##To make a cumulative histogram to show unique annotations with respect to crude and previous fractions

# Define the crude sample (or samples) to act as the baseline                     ###USER-INPUT###
remove_sample <- "HGMA_1598" 

# Filter out rows containing the specified sample
#Annotations.with.samples <- Annotations.with.samples #%>%                        ## Remove #'s if wanting to filter out crude or anotehr sample
#  filter(!str_detect(Samples, remove_sample))                                    ## Same here

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

ggsave(filename = paste0(folder, "/cumulative-histogram.svg"), 
       plot = plot_obj3, 
       width = 50, 
       height = 40, 
       units = "cm")

                #BUBBLE CHART: FOR MS2QUERY                                
                library(ggplot2)
                library(reshape2) 
                library(dplyr)
                library(tidyr)
                library(readr)
                library(svglite)
                library(readxl)
                
                ##---##
        
                metadata <- read.csv("HGM/A - Analysis.csv")
                names(metadata)[1] <- "Sample"
                metadata <- select(metadata, Sample, Book.Code)
                
                #Read in CSV's and merge
                ms2query.data <- paste0(folder, "/ms2query/ms2query.csv")
                class.data <- read.csv(ms2query.data) %>%                            #####!!!###### Check spelling of filename
                  filter(ms2query_model_prediction >= 0.63)
                names(class.data)[10] <- "row.ID"
                class.data <- select(class.data, "row.ID", "npc_superclass_results") %>%
                  filter(npc_superclass_results != "None")
                
                quant.data <- paste0(folder, "/mzmine/DATA_iimn_gnps_quant.csv")
                quant.data <- read.csv(quant.data)                    #####!!!######  Check spelling of filename
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
                
                ggsave(filename = paste0(folder, "/bubblechart.svg"), 
                       plot = plot_obj, 
                       width = 50, 
                       height = 40, 
                       units = "cm")
