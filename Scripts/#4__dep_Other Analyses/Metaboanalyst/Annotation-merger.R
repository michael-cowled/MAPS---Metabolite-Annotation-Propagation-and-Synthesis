## Spectral library matches to Quant File for metaboanalyst

library(dplyr)
library(tidyr)
library(stringr)

#Import Quant and ms1+ms2 csv's:
quant.data <- read.csv("DATA_iimn_gnps_quant.csv")                    #####!!!######  Check spelling of filename
quant.data <- quant.data[, -c(2:13)] # Removes unnecessary columns
colnames(quant.data)[1] <- "feature.ID"

annotation.data <- read.csv("final-annotation-df.csv") 
annotation.data <- select(annotation.data, feature.ID, Best.Annotation) %>%
  filter(!is.na(Best.Annotation))

metadata <- read.csv("metadata.csv")

merged.data <- annotation.data %>%
  left_join(quant.data, by = "feature.ID") %>%
  select(-X, -feature.ID) %>%
  t()
  metaboanalyst.data <- cbind(filename = row.names(merged.data), merged.data)
  metaboanalyst.data <- as.data.frame(metaboanalyst.data)
  metaboanalyst.data$filename <- str_remove(metaboanalyst.data$filename, "\\.mzML\\.Peak\\.area$")

  metaboanalyst.data <-  metaboanalyst.data %>%
    left_join(metadata, by = "filename") ##Append metadata
  
  ##Reorder columns such that metadata is second:
  col_names <- colnames(metaboanalyst.data)
  new_order <- c(col_names[1], col_names[length(col_names)], col_names[2:(length(col_names) - 1)])
  metaboanalyst.data <- metaboanalyst.data[, new_order]
  metaboanalyst.data[1,2] <- "Metadata"
  
  #Create file to import into metaboanalyst

  write.csv(metaboanalyst.data, "metaboanalyst-tidied.csv")      ## Note, just remove first column and row before importing