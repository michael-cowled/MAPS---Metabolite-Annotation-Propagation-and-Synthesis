library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)

#Importing data
#note data file is already pre-tidied with the relatively-empty columns deleted before importing
#m/z and retention time copied to a new table titled "feature_id_to_mass" with the corresponding "ID"

data <- read.csv("C:/Users/miccow/Desktop/Apple_Fungi_Metabolic_Profiling/MZMINE full dataset/GNPS_tidied.csv")
metadata <- read.csv("C:/Users/miccow/Desktop/Apple_Fungi_Metabolic_Profiling/MZMINE full dataset/metadata.csv")

#Remove primary metabolites from total dataset (data) by using tidied app.feats df as list of viable features
colnames(data)[1] <- "feature.ID"
simple.app <- app.feats[,c(1,5)]
data2 <- merge(simple.app, data,
               by = 'feature.ID', all = TRUE)
data2 <- filter(data2, !is.na(componentindex))
data2 <- select(data2, -componentindex)

#Tidying data
col.names <- colnames(data2) # makes a list of colnames
col.names <- gsub("\\..*","",col.names) # tidies colnames
colnames(data2) <- col.names # updates colnames in data
data2 <- t(data2) # transposes data

## Generates feature.ID.2 from data2
data3 <- filter(as.data.frame(data2), V1==1.0000)
data3 <- t(data3) # transposes data back
data3 <- as.data.frame(data3)
row.ID2 <- rownames(data3) ## note, a different row.ID to the other one used
data3 <- cbind(row.ID2, data3)
colnames(data3) <- c("row.ID2", "feature.ID")

#Tidying metadata
metadata2 <- separate(metadata, filename, c("filename", "discard"))
rownames(metadata2) <- metadata2[,1]

#Merging with the metadata
data.merged <- merge(metadata2, data2,
                     by = 'row.names', all = TRUE) %>%
  select(-Row.names, -discard) %>%
  filter(filename != "PDAd7con")  ## PDAd7con removed due to being contaminated ##

##Blanks to be filtered out.
blanks <- filter(data.merged, data.merged$ATTRIBUTE_Type == "Blank_Type")
blank.sums <- summarise_all(blanks, ~if(is.numeric(.)) sum(.) else "Total")
n <- ncol(blank.sums)
col.to.be.filtered <- 8
data.filt.blanks <- data.merged
list.of.blank.features <- list()

##Finds and filters blank features
while (col.to.be.filtered <= n) {
  if (blank.sums[col.to.be.filtered] >= 1000) { ##Threshold of 1000 intensity across all blanks
    list.of.blank.features <- append(list.of.blank.features, col.to.be.filtered)
  } 
  col.to.be.filtered <- col.to.be.filtered + 1
  print(col.to.be.filtered)
}

data.filt.blanks <- data.filt.blanks[,-c(unlist(list.of.blank.features))]
data.filt.blanks <- filter(data.filt.blanks, ATTRIBUTE_Type != "Blank_Type") %>%
  filter(V2 >= 0) ## removes blank datasets

paste0(length(list.of.blank.features), " blank features found and removed!")


##-------------------------------------##

## Splitting into different categories and finding which features are present
apple <- filter(data.filt.blanks, ATTRIBUTE_Type == "Apple_natural" | ATTRIBUTE_Type == "Apple_excision") #Can edit if only interested in excisions of whole black apples
apple.feat.only <- apple[,-c(1:7)]
apple.sums <- colSums(apple.feat.only)
lab <- filter(data.filt.blanks, ATTRIBUTE_Type == "Laboratory")
lab.feat.only <- lab[,-c(1:7)]
lab.sums <- colSums(lab.feat.only)

paste0(length(data.filt.blanks), " features in both apple and lab")
paste0(sum(apple.sums > 0), " features in apple")
paste0(sum(lab.sums > 0), " features in lab")
paste0(sum(lab.sums > 0 & apple.sums > 0), " features in both apple and lab")

## Correlation to featureIDs
featureIDs <- read.csv("C:/Users/miccow/Desktop/Apple_Fungi_Metabolic_Profiling/MZMINE full dataset/feature_id_to_mass_pre-tidied.csv")
in.both <- lab.sums > 0 & apple.sums > 0
in.both <- as.data.frame(in.both)
row.ID2 <- rownames(in.both) ##not comparable with the above data3 object
in.both <- cbind(row.ID2,in.both)
in.both <- merge(data3, in.both, by = "row.ID2")
in.both <- select(in.both, -row.ID2) # removes temporary row.ID
merged.df <- merge(featureIDs, in.both, by = "feature.ID") ## now feature.ID is compatible.

# Now some analysis of features in.both
merged.df.with.metadata2 <- merge(app.feats, merged.df, by = "feature.ID")
in.both2 <- filter(merged.df.with.metadata2, in.both == TRUE) 

# Piechart of NPclassifier superclasses

df <- in.both2 %>% 
  group_by(NPC_NPC.superclass) %>% # Variable to be transformed
  count() %>% 
  ungroup() %>% 
  filter(!is.na(NPC_NPC.superclass)) %>%
  mutate(perc = `n` / sum(`n`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(round(perc, 3))) %>%
  mutate(lab.ypos = cumsum(perc) - 0.5*perc)
df$x <- paste0(df$NPC_NPC.superclass, " (", df$labels, ")")

df2 <- df %>% 
  mutate(csum = rev(cumsum(rev(perc*100))), 
         pos = (perc*100)/2 + lead(csum, 1),
         pos = if_else(is.na(pos), (perc*100)/2, pos))

no_of_colors <- nrow(df)
mycols <- viridis_pal(option = "A")(no_of_colors)

ggplot(df, aes(x = "", y = perc, fill = x)) +
  geom_bar(width = 3, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  scale_fill_manual(values = mycols) +
  theme_void() +
  guides(fill=guide_legend(title="NP Superclass"))