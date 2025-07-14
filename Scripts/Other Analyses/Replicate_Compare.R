### META_COMPARE.R ###
##To compare metadata

#metadata - GNPS compatible file (tab delimited .txt file)
metadata <- read_tsv("metadata.txt")
metadata$filename <- sub("\\.mzML$", "", metadata$filename)
metadata <- as.data.frame(metadata)
colnames(metadata)[2] <- "metadata"

R1 <- filter(metadata, Comparison == "R1")
R2 <- filter(metadata, Comparison == "R2")
R3 <- filter(metadata, Comparison == "R3")

final.annotation.df3 <- final.annotation.df2
final.annotation.df4 <- filter(final.annotation.df2, !is.na(Best.Annotation))

# Initialize the new column with FALSE
final.annotation.df4$ContainsR1Filename <- FALSE

# Loop through each filename in R1$filename
for (i in 1:nrow(R1)) {
  filename <- R1$filename[i] 
  print(filename)
  final.annotation.df4$ContainsR1Filename <- final.annotation.df4$ContainsR1Filename | grepl(filename, final.annotation.df4$Samples)
}

# Initialize the new column with FALSE
final.annotation.df4$ContainsR2Filename <- FALSE

# Loop through each filename in R1$filename
for (i in 1:nrow(R2)) {
  filename <- R2$filename[i] 
  print(filename)
  final.annotation.df4$ContainsR2Filename <- final.annotation.df4$ContainsR2Filename | grepl(filename, final.annotation.df4$Samples)
}

# Initialize the new column with FALSE
final.annotation.df4$ContainsR3Filename <- FALSE

# Loop through each filename in R1$filename
for (i in 1:nrow(R3)) {
  filename <- R3$filename[i] 
  print(filename)
  final.annotation.df4$ContainsR3Filename <- final.annotation.df4$ContainsR3Filename | grepl(filename, final.annotation.df4$Samples)
}

#Tidy up for features that are relevant for the venn digaram 
final.annotation.df4 <- final.annotation.df4 %>%
  mutate(ContainsAny = ContainsR3Filename + ContainsR2Filename + ContainsR1Filename) %>%
  filter(ContainsAny >=1)


library(tidyverse)
#install.packages("ggvenn")
library(ggvenn)
#> Loading required package: grid
Name <- c("R1", "R2", "R3")

final.annotation.df4 %>%
  ggplot() +
  geom_venn(aes(A = ContainsR1Filename, B = ContainsR2Filename, C = ContainsR3Filename),
            set_names = Name, fill_color = c("black", "black", "black"), text_color = "white") +
  theme_void()