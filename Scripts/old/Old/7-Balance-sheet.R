##Balance Sheet of Annotations
#How to extract useful metrics for each fraction (number of annotations, and of each category)

#Sample based:
#Read in either the most recently generated or that manually copied to documents:
final.annotation.df <- read.csv("final-annotation-df.csv")  #total number of features

##Any filtering for specific samples?
filter_sample <- "HGMA_1623"

# Filter out rows containing the specified sample
final.annotation.df <- final.annotation.df %>%
  filter(!str_detect(Samples, filter_sample ))

num_lv1 <- filter(final.annotation.df, !is.na(authentic.standard))
num_lv2 <- filter(final.annotation.df, is.na(authentic.standard), !is.na(Best.Annotation.with.Analogues))
num_lv3 <- filter(final.annotation.df, 
                  is.na(Best.Annotation.with.Analogues) & 
                    (!is.na(canopus.NPC.superclass) | !is.na(ms2query.NPC.superclass)))

nrow(final.annotation.df)
nrow(num_lv1)
nrow(num_lv2)
nrow(num_lv3)
nrow(final.annotation.df)-nrow(num_lv3)-nrow(num_lv2)-nrow(num_lv1)





##Ongoing problem for another time------
#Combine across two different datasets to show which are obtained from which analysis
#now we have two datasets: the AQ run in Phe-Hex and DCM in Phe-Hex) – 
#how to combine? What identifier as names are not working? Very important for when we try to include other analysis platforms

#metdata based:
metadata <- read.csv("metadata.csv")

#Read in either the most recently generated or that manually copied to documents:
final.annotation.df <- read.csv("final-annotation-df.csv")  #total number of features

#Specific Metadata filtering
Sample.filter <- filter(metadata, Group == "MC01-61" & ID2 == "Female")

# Initialize the new column with FALSE
balance.sheet <- final.annotation.df
balance.sheet$ContainsR1Filename <- FALSE

# Loop through each filename in R1$filename
for (i in 1:nrow(Sample.filter)) {
  filename <- Sample.filter$filename[i] 
  print(filename)
  balance.sheet$ContainsR1Filename <- balance.sheet$ContainsR1Filename | grepl(filename, balance.sheet$Samples)
}

nrow(balance.sheet)
nrow(num_lv1)
nrow(num_lv2)
nrow(num_lv3)
nrow(balance.sheet)-nrow(num_lv3)-nrow(num_lv2)-nrow(num_lv1)