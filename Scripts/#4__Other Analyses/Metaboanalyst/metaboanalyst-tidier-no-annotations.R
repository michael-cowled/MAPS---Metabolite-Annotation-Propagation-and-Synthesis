## Spectral library matches to Quant File for metaboanalyst

library(dplyr)
library(tidyr)
library(stringr)

#Import Quant:

quant.data <- read.csv("DATA_iimn_gnps_quant.csv")                    #####!!!######  Check spelling of filename
quant.data <- quant.data[, -c(2:13)] # Removes unnecessary columns

cfucount.data <- read.csv("cfucount.csv")

metadata <- read.csv("metadata.csv")

merged.data <- quant.data %>%                   ##Actually no merging occurs, but kept the same as Annotation-merger.R
  mutate(across(3:ncol(.), ~ ifelse(. < 1000, 0, .)))  #Applies a filter to remove features potentially below the noise

merged.data <- merged.data[!is.na(merged.data[, 3]), ]

#Subtract blank values
#First: Avg of Media blanks and instrument blanks
#merged.data$Media_avg <- rowMeans(merged.data[, c("HGMA_0461.mzML.Peak.area", "HGMA_0462.mzML.Peak.area", "HGMA_0463.mzML.Peak.area")], na.rm = TRUE)
merged.data$Blank_avg <- rowMeans(merged.data[, c("Prerun_blank_01.mzML.Peak.area", "Prerun_blank_02.mzML.Peak.area", "Prerun_blank_03.mzML.Peak.area", "Prerun_blank_04.mzML.Peak.area", "Prerun_blank_05.mzML.Peak.area",
                                                  "Postrun_blank_01.mzML.Peak.area", "Postrun_blank_02.mzML.Peak.area", "Postrun_blank_03.mzML.Peak.area", "Postrun_blank_04.mzML.Peak.area", "Postrun_blank_05.mzML.Peak.area")], na.rm = TRUE)

#Second: subtract instrument blanks
merged.data.subtract <- merged.data
merged.data.subtract[, 3:ncol(merged.data)] <- merged.data[, 3:ncol(merged.data)] - merged.data$Blank_avg

#Third: subtract media blanks, making sure not to double subtract an instrument blank
##merged.data.subtract$Media_avg[merged.data.subtract$Media_avg < 0] <- 0
##merged.data.subtract[, 3:ncol(merged.data)] <- merged.data[, 3:ncol(merged.data)] - merged.data$Media_avg

#Fourth: remove columns associated with media and blanks
merged.data.subtract <- select(merged.data.subtract, -"Blank_avg", 
                               -"Prerun_blank_01.mzML.Peak.area", -"Prerun_blank_02.mzML.Peak.area", -"Prerun_blank_03.mzML.Peak.area", -"Prerun_blank_04.mzML.Peak.area", -"Prerun_blank_05.mzML.Peak.area",
                               -"Postrun_blank_01.mzML.Peak.area", -"Postrun_blank_02.mzML.Peak.area", -"Postrun_blank_03.mzML.Peak.area", -"Postrun_blank_04.mzML.Peak.area", -"Postrun_blank_05.mzML.Peak.area", 
                               -"Midrun_blank_01.mzML.Peak.area", -"Midrun_blank_02.mzML.Peak.area", -"HGMA_0635.mzML.Peak.area", -"HGMA_0654.mzML.Peak.area", -"HGMA_0764_A.mzML.Peak.area", -"HGMA_0764_B.mzML.Peak.area", -X)

#Lastly: Correct all negative values to be 0.
merged.data.subtract <- merged.data.subtract %>%
  mutate(across(3:ncol(.), ~ ifelse(. < 0, 0, .)))

##The following applies to filtering out features with no annotation to metaboanalyst
metaboanalyst.data <- merged.data.subtract[!is.na(merged.data.subtract[, 2]), ] %>%
  select(-row.ID) %>%
  t()
metaboanalyst.data <- cbind(filename = row.names(metaboanalyst.data), metaboanalyst.data)
metaboanalyst.data <- as.data.frame(metaboanalyst.data)
metaboanalyst.data$filename <- str_remove(metaboanalyst.data$filename, "\\.mzML\\.Peak\\.area$")

#Find and remove features which are zero across all samples
metaboanalyst.data[, 3:ncol(metaboanalyst.data)] <- lapply(metaboanalyst.data[, 3:ncol(metaboanalyst.data)], as.numeric)
# Calculate the column sums for columns 3 onwards
col_sums <- colSums(metaboanalyst.data[, 3:ncol(metaboanalyst.data)], na.rm = TRUE)
# Identify columns with sums greater than zero
columns_to_keep <- col_sums > 0
# Keep the first two columns and columns with sums greater than zero
metaboanalyst.data <- metaboanalyst.data[, c(TRUE, TRUE, columns_to_keep)]

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
export.file <- filter(export.file, filename != export.file[1,1])
names(export.file)[1:2] <- c("filename", "Metadata")

#Create file to import into metaboanalyst

write.csv(export.file, "metaboanalyst-tidied.csv")      ## Note, just remove first column before importing