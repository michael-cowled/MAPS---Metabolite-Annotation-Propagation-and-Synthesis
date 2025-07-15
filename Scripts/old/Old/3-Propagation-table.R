###Propagation Table
##Using full.annotation.data, perform propagations

##Create a summary.annotation.data table, that is filtered for all annotations satisfying a certain degree of confidence
summary.annotation.data <- full.annotation.data
#Apply the conditions
summary.annotation.data$csi.compound.name[summary.annotation.data$csi.confidence.score < csi.prob] <- NA
summary.annotation.data$csi.smiles[summary.annotation.data$csi.confidence.score < csi.prob] <- NA
summary.annotation.data$ms2query.analogue.compound.name[summary.annotation.data$ms2query.score < ms2query.prob] <- NA
summary.annotation.data$ms2query.analogue.compound.name[summary.annotation.data$ms2query.mzdiff > 0.001] <- NA

# Initialize columns in propagated.annotation.data
propagated.annotation.data <- full.annotation.data %>%
  mutate(
    Probable.Analogue.Of = NA,
    Propagated.Feature.ID = NA,
    Propagated.Annotation.Type = NA  # Initialize the new column
  )

# Function to get the first non-NA value from specified columns and the column name
get_result <- function(paired_value, summary_data) {
  # Define the columns to check for non-NA values
  columns_to_check <- c("gnps.compound.name", "csi.compound.name", 
                        "ms2query.analogue.compound.name")
  
  # Filter the relevant row and select the columns
  data_subset <- summary_data %>%
    filter(feature.ID == paired_value) %>%
    select(all_of(columns_to_check), csi.smiles)
  
  # Find the first non-NA and non-'null' value
  result <- NA
  column_name <- NA
  for (col in columns_to_check) {
    value <- data_subset[[col]]
    if (!is.na(value) && value != 'null') {
      result <- value
      column_name <- col
      break
    }
  }
  
  # If no valid result is found, fallback to 'csi.smiles'
  if (is.na(result) || result == 'null') {
    csi_value <- data_subset$csi.smiles
    if (!is.na(csi_value) && csi_value != 'null') {
      result <- csi_value
      column_name <- "csi.smiles"
    }
  }
  
  # Return the result and corresponding column name
  list(value = result, column = column_name)
}

# Function: Finds the features linked by an edge in a GNPS cluster
paired.feature.finder <- function(ID) {
  filtered_pairs <- gnps.cluster.pairs %>%
    filter(CLUSTERID1 == ID | CLUSTERID2 == ID) %>%
    arrange(desc(Cosine))  # Sort by decreasing Cosine
  
  # Extract paired values
  paired_values <- filtered_pairs %>%
    mutate(paired_value = ifelse(CLUSTERID1 == ID, CLUSTERID2, CLUSTERID1)) %>%
    pull(paired_value)
  
  return(paired_values)
}

# Identify the rows where all specified columns contain NA
na.rows <- summary.annotation.data %>%
  filter(
    is.na(gnps.compound.name) &
      is.na(csi.compound.name) &
      is.na(csi.smiles) &
      is.na(ms2query.analogue.compound.name)
  )

# Extract the 'feature.ID' column from these rows
na.feature.ids <- na.rows$feature.ID

# Iterate over na.feature.ids to populate propagated.annotation.data
for (i in na.feature.ids) {
  paired_values <- paired.feature.finder(i)
  print(paste("Paired values for feature ID", i, ":", paste(paired_values, collapse = ", ")))
  
  # Initialize variables to store results and corresponding paired value
  result_data <- list(value = NA, column = NA)
  selected_paired_value <- NA
  
  for (value in paired_values) {
    result_data <- get_result(value, summary.annotation.data)
    if (!is.na(result_data$value)) {
      selected_paired_value <- value  # Store the paired value used for the result
      break  # Exit the loop if a valid result is found
    }
  }
  
  # Update the new columns in propagated.annotation.data
  propagated.annotation.data <- propagated.annotation.data %>%
    mutate(
      Probable.Analogue.Of = case_when(
        feature.ID == i ~ result_data$value,
        TRUE ~ Probable.Analogue.Of
      ),
      Propagated.Feature.ID = case_when(
        feature.ID == i ~ selected_paired_value,
        TRUE ~ Propagated.Feature.ID
      ),
      Propagated.Annotation.Type = case_when(
        feature.ID == i ~ result_data$column,
        TRUE ~ Propagated.Annotation.Type
      )
    )
}

#Calculating the mz difference of the analogue to the propagated feature:
# Convert mz column to numeric if it's not already
if (!is.numeric(summary.annotation.data$mz)) {
  summary.annotation.data$mz <- as.numeric(summary.annotation.data$mz)
}

# Add a new column for the m/z difference
propagated.annotation.data <- propagated.annotation.data %>%
  mutate(Propagated.Annotation.mz.Diff = NA)

# Iterate over the rows of propagated.annotation.data
for (i in 1:nrow(propagated.annotation.data)) {
  # Get the feature ID and the propagated feature ID
  feature_id <- propagated.annotation.data$feature.ID[i]
  propagated_feature_id <- propagated.annotation.data$Propagated.Feature.ID[i]
  
  # If a propagated feature ID exists
  if (!is.na(propagated_feature_id)) {
    # Get the m/z values for both feature IDs
    mz_value <- summary.annotation.data$mz[summary.annotation.data$feature.ID == feature_id]
    propagated_mz_value <- summary.annotation.data$mz[summary.annotation.data$feature.ID == propagated_feature_id]
    
    # Calculate the m/z difference
    mz_diff <- mz_value - propagated_mz_value
    
    # Update the new column
    propagated.annotation.data$Propagated.Annotation.mz.Diff[i] <- mz_diff
  }
}

# Inside the loop:
if (!is.na(propagated_feature_id)) {
  mz_value <- summary.annotation.data$mz[summary.annotation.data$feature.ID == feature_id]
  propagated_mz_value <- summary.annotation.data$mz[summary.annotation.data$feature.ID == propagated_feature_id]
  
  if (!is.na(mz_value) && !is.na(propagated_mz_value)) {  # Check for NA
    mz_diff <- mz_value - propagated_mz_value
    propagated.annotation.data$Propagated.Annotation.mz.Diff[i] <- mz_diff
  }
}

##Now to append sample.list corresponding to each feature

# Select columns containing ".area" and rename the first column
sample.data <- sample.data %>%
  select(id, contains(".area"))

# Rename the first column to "feature.ID"
colnames(sample.data)[1] <- "feature.ID"

# Get the column names of the data frame
colnames_sample <- colnames(sample.data)

# Remove the prefix "datafile." and the suffix ".mzML.Peak.area" if present
colnames_sample <- sub("^datafile\\.", "", colnames_sample)  # Remove the prefix
colnames_sample <- sub("\\.mzML\\.area$", "", colnames_sample)  # Remove the suffix

# Replace all positive areas with 1 i.e. presence.absence
colnames(sample.data) <- colnames_sample
sample.data[, 2:ncol(sample.data)] <- lapply(sample.data[, 2:ncol(sample.data)], function(x) {
  x[x > 0] <- 1
  return(x)
})

# Create a new column with presence list separated by semicolons
sample.data$Samples <- apply(sample.data[, 2:ncol(sample.data)], 1, function(row) {
  paste(colnames(sample.data[, 2:ncol(sample.data)])[which(row == 1)], collapse = "; ")
})

# Select only the first and the last columns (which includes the new 'Samples' column)
sample.data <- sample.data[, c(1, ncol(sample.data))]

#Append Sample data
propagated.annotation.data.with.samples <- propagated.annotation.data %>%
  full_join(sample.data, by = "feature.ID") 

##Collapsing Ion Identity Networks - performed in such a way to retain the best annotation (if multiple)

#A new editing df where features are sequentially removed:
propagation.df <- propagated.annotation.data.with.samples

#A new df with "good" annotations extracted in a similar manner to unpropagated version.
new.summary.annotation.data <- propagated.annotation.data.with.samples
#Apply the conditions
new.summary.annotation.data$csi.compound.name[new.summary.annotation.data$csi.confidence.score < csi.prob] <- NA
new.summary.annotation.data$csi.smiles[new.summary.annotation.data$csi.confidence.score < csi.prob] <- NA
new.summary.annotation.data$ms2query.analogue.compound.name[new.summary.annotation.data$ms2query.score < ms2query.prob] <- NA

iin.features <- filter(new.summary.annotation.data, !is.na(ion.identity.ID))

# Define a function to process each IIN group
process_iin_group <- function(group_data, feature_columns) {
  # Initialize a vector to track which rows to retain
  rows_to_retain <- logical(nrow(group_data))
  
  for (col in feature_columns) {
    if (any(!is.na(group_data[[col]]))) {
      values <- group_data$feature.ID[!is.na(group_data[[col]])]
      if (length(values) > 1) {
        # Multiple values meet the condition, retain the lowest values
        highest_value <- max(values, na.rm = TRUE)
        rows_to_retain <- group_data$feature.ID %in% values & group_data$feature.ID != highest_value
      } else {
        # Only one value meets the condition, retain this row
        rows_to_retain <- group_data$feature.ID == values
      }
      break  # Exit loop once a condition is met
    }
  }
  
  if (all(!rows_to_retain)) {
    # No rows were marked to retain, so retain the one with the minimum feature.ID for fallback
    rows_to_retain <- group_data$feature.ID == min(group_data$feature.ID, na.rm = TRUE)
  }
  
  return(group_data[rows_to_retain, ])
}

# Define a function to update the main data frame
update_data_frame <- function(df, results, id_column) {
  df_non_na <- df %>% filter(!is.na(!!sym(id_column)))
  results_df <- bind_rows(results)
  
  updated_df <- df_non_na %>% filter(feature.ID %in% results_df$feature.ID)
  df_na <- df %>% filter(is.na(!!sym(id_column)))
  
  final_df <- bind_rows(updated_df, df_na)
  return(final_df)
}

# Main processing function
process_all_features <- function(data, id_column, feature_columns) {
  unique_ids <- unique(data[[id_column]][!is.na(data[[id_column]])])
  results <- lapply(unique_ids, function(id) {
    group_data <- filter(data, !!sym(id_column) == id)
    process_iin_group(group_data, feature_columns)
  })
  names(results) <- as.character(unique_ids)
  return(results)
}

# Print results for each ID
for (id in names(results)) {
  cat("Results for ion.identity.ID =", id, ":\n")
  print(results[[id]])
}

# Update propagation.df by removing rows based on results and retaining NA IDs
propagation.df <- update_data_frame(propagation.df, results, "ion.identity.ID")