###Meta-Compare###
#Aims to compare datasets or samples by metadata

library(ggplot2)
library(reshape2) 
library(dplyr)
library(tidyr)
library(readr)

#GNPS2 Output Link (Private)

gnps.link <- "https://gnps2.org/status?task=e24e2a851f5d4aa79e4c4f93fa5603f7" 

# GNPS2 - works for v0.1.2, no metadata required
gnps.link <- gsub("status", "resultfile", gnps.link) # Updates the user-input link to direct to "resultfile"
gnps.annotation.data <- read_tsv(paste0(gnps.link, "&file=nf_output/library/merged_results_with_gnps.tsv"))
gnps.annotation.data <- gnps.annotation.data[, c(2, 15)]
names(gnps.annotation.data) <- c("feature.ID", "gnps.compound.name")

# MZmine
quant.data <- read.csv("DATA_iimn_gnps_quant.csv")                    #####!!!######  Check spelling of filename
quant.data <- quant.data[, -c(2:13)] # Removes unneccessary columns
colnames_quant <- colnames(quant.data) # Get the column names of quant.data
colnames_quant <- sub("\\.mzML\\.Peak\\.area$", "", colnames_quant) # Remove the suffix ".mzML.Peak.area" if present
colnames(quant.data) <- colnames_quant # Assign the modified column names back to quant.data
names(quant.data)[1] <- 'feature.ID'

#metadata - GNPS compatible file (tab delimited .txt file)
metadata <- read_tsv("metadata.txt")
metadata$filename <- sub("\\.mzML$", "", metadata$filename)
metadata <- as.data.frame(metadata)
colnames(metadata)[2] <- "metadata"

#Merge GNPS and MZMINE
merged.data <- right_join(quant.data, gnps.annotation.data, by = "feature.ID") ## Reduces to the number of annotated featured only
merged.data <- t(merged.data)
filenames <- rownames(merged.data) ## note, a different row.ID to the other one used
merged.data <- cbind(filenames, merged.data)
colnames(merged.data)[1] <- "filename"
merged.data <- as.data.frame(merged.data)

#Append metadata
merged.data.with.metadata <- full_join(metadata, merged.data, by = "filename")

#For metadata_attribute_1:
m.a.1 <- filter(merged.data.with.metadata, metadata == "MC01-55") %>%
  select(-filename, -metadata)
m.a.1[] <- lapply(m.a.1, as.numeric)
m.a.1 <- colSums(m.a.1, na.rm = TRUE)
m.a.1 <- as.data.frame(m.a.1)
arb <- rownames(m.a.1)
arb <- as.data.frame(arb)
m.a.1 <- cbind(arb, m.a.1)

#For metadata_attribute_2:
m.a.2 <- filter(merged.data.with.metadata, metadata == "MC01-57") %>%
  select(-filename, -metadata)
m.a.2[] <- lapply(m.a.2, as.numeric)
m.a.2 <- colSums(m.a.2, na.rm = TRUE)
m.a.2 <- as.data.frame(m.a.2)
arb <- rownames(m.a.2)
arb <- as.data.frame(arb)
m.a.2 <- cbind(arb, m.a.2)

##Now to compare:
annotations <- merged.data[nrow(merged.data),] %>%
  select(-filename)
annotations <- t(annotations)
arb <- rownames(annotations)
annotations <- cbind(arb, annotations)
annotations <- as.data.frame(annotations)

#merge datasets
merged.annotations <- right_join(annotations, m.a.1, by = "arb")
merged.annotations <- right_join(merged.annotations, m.a.2, by = "arb")

write.csv(merged.annotations, "metadata-compared.csv")