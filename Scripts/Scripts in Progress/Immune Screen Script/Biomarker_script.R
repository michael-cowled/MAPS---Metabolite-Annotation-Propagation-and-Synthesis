###BIOMARKER SCRIPT###
##To use Annotations, Samples and Immune Screen Results combined.

##First: To Combine Annotations and Samples
final.annotation.df3 <- final.annotation.df2 %>%
  select(feature.ID, mz, Best.Annotation.with.Analogues, Best.Annotation.with.Analogues.Smiles, confidence.level)

#Append to Samples
samples.df.with.annotations <- samples.df %>%
  full_join(final.annotation.df3, by = "feature.ID")

##Second: Optional df to extract top 10 best features representative of a sample.
samples.df.with.annotations 

top_10_features <- samples.df.with.annotations %>%
  group_by(samples) %>%
  top_n(10, area) %>%
  ungroup()

head(top_10_features) # Display the top of the resulting data frame

###Incoporated in Annotations Script