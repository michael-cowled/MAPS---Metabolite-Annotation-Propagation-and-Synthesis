## Spectral library matches to Quant File for metaboanalyst

library(dplyr)
library(tidyr)
library(stringr)

#Import Quant and ms1+ms2 csv's:
annotation.data <- read.csv("ms1-and-ms2.csv")                    #####!!!###### Check spelling of filename
annotation.data <- annotation.data[, c(1, 32)] %>%
  mutate(across(everything(), as.character)) %>%
  mutate(across(everything(), ~ na_if(., "")))
names(annotation.data) <- c("row.ID", "Annotation")
annotation.data$row.ID <- as.numeric(annotation.data$row.ID)

quant.data <- read.csv("DATA_iimn_gnps_quant.csv")                    #####!!!######  Check spelling of filename
quant.data <- quant.data[, -c(2:13)] # Removes unnecessary columns

cfucount.data <- read.csv("cfucount.csv")

metadata <- read.csv("metadata.csv")

merged.data <- annotation.data %>%
  left_join(quant.data, by = "row.ID") %>%
  select(-X)

merged.data <- merged.data[!is.na(merged.data[, 3]), ]

##The following applies to filtering out features with no annotation to metaboanalyst
metaboanalyst.data <- merged.data.subtract[!is.na(merged.data.subtract[, 2]), ] %>%
  select(-row.ID) %>%
  t()
  metaboanalyst.data <- cbind(filename = row.names(metaboanalyst.data), metaboanalyst.data)
  metaboanalyst.data <- as.data.frame(metaboanalyst.data)
  metaboanalyst.data$filename <- str_remove(metaboanalyst.data$filename, "\\.mzML\\.Peak\\.area$")

  
#Normalise to CFU
  #First: verify filenames in quant are the same as those in cfucount
  metaboanalyst.data2 <- metaboanalyst.data %>%
    left_join(cfucount.data, by = "filename") ##Append CFUdata
  Annotations <- metaboanalyst.data2[1,] ##Records annotations
  metaboanalyst.data2 <- metaboanalyst.data2[-1, ]  ##Temporarily removes annotations for ease of manipulation
  
  #Second: Normalise by dividing by cfu
  metaboanalyst.data2[, -1] <- lapply(metaboanalyst.data2[, -1], as.numeric)
  
  # Divide all columns (except row 1 and column 1) by metaboanalyst.data2$cfumL
  metaboanalyst.data2[, 2:ncol(metaboanalyst.data2)] <- metaboanalyst.data2[, 2:ncol(metaboanalyst.data2)] / metaboanalyst.data2$CFUml

  metaboanalyst.data2_filtered <- metaboanalyst.data2 %>%
    filter(!is.na(CFUml))
  export.file <- rbind(Annotations,metaboanalyst.data2_filtered) %>%
    select(-CFUml) %>%
    left_join(metadata, by = "filename") ##Append metadata
  
  ##Reorder columns such that metadata is second:
  col_names <- colnames(export.file)
  new_order <- c(col_names[1], col_names[length(col_names)], col_names[2:(length(col_names) - 1)])
  export.file <- export.file[, new_order]
  export.file[1,2] <- "Metadata"
  
  #Create file to import into metaboanalyst

  write.csv(export.file, "metaboanalyst-tidied.csv")      ## Note, just remove first column and row before importing