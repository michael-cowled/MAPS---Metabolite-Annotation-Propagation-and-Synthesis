####Cytoscape Tidier####

##Uses output from Annotation-Table.R##

cytoscape <- read.csv("cytoscape.csv")

level1.annotation <- filter(propagation.df, !is.na(authentic.standard)) %>%
  select(feature.ID, authentic.standard, authentic.standard.smiles)
names(level1.annotation) <- c("feature.ID", "annotation", "smiles")
level1.annotation$confidence.level <- "1"
not.level1.annotation <- filter(propagation.df, !feature.ID %in% level1.annotation$feature.ID)

gnps <- filter(not.level1.annotation, !is.na(gnps.compound.name)) %>%
  select(feature.ID, gnps.compound.name, gnps.smiles)
names(gnps) <- c("feature.ID", "annotation", "smiles")
gnps$confidence.level <- "2"
not.gnps <- filter(not.level1.annotation, !feature.ID %in% gnps$feature.ID)

csi <- filter(not.gnps, !is.na(csi.compound.name) | csi.compound.name != "null") %>%
  filter(csi.confidence.score > csi.prob) %>%
  select(feature.ID, csi.compound.name, csi.smiles)
names(csi) <- c("feature.ID", "annotation", "smiles")
csi$confidence.level <- "2"
not.csi <- filter(not.gnps, !feature.ID %in% csi$feature.ID)

ms2q.real <- filter(not.csi, !is.na(ms2query.analogue.compound.name)) %>%
  filter(ms2query.score > ms2query.prob) %>%
  filter(ms2query.mzdiff <= 0.001) %>%
  select(feature.ID, ms2query.analogue.compound.name, ms2query.smiles)
names(ms2q.real) <- c("feature.ID", "annotation", "smiles")
ms2q.real$confidence.level <- "2"
not.real.ms2q <- filter(not.csi, !feature.ID %in% ms2q.real$feature.ID)

ms2q <- filter(not.real.ms2q, !is.na(ms2query.analogue.compound.name)) %>%
                 filter(ms2query.score > ms2query.prob) %>%
  select(feature.ID, ms2query.analogue.compound.name, ms2query.smiles)
names(ms2q) <- c("feature.ID", "annotation", "smiles")
ms2q$annotation <- paste0("Analogue of ", ms2q$annotation)
ms2q$confidence.level <- "3"
not.ms2q <- filter(not.real.ms2q, !feature.ID %in% ms2q$feature.ID)

cytoscape.annotations <- rbind(level1.annotation, gnps)
cytoscape.annotations <- rbind(cytoscape.annotations, csi)
cytoscape.annotations <- rbind(cytoscape.annotations, ms2q.real)

final.annotation.df <-  left_join(propagation.df, cytoscape.annotations, by = "feature.ID") %>%
  select(-confidence.level)
names(final.annotation.df)[names(final.annotation.df) == "annotation"] <- "Best.Annotation"
names(final.annotation.df)[names(final.annotation.df) == "smiles"] <- "Best.Annotation.Smiles"

cytoscape.annotations <- rbind(cytoscape.annotations, ms2q)

final.annotation.df <-  left_join(final.annotation.df, cytoscape.annotations, by = "feature.ID")
names(final.annotation.df)[names(final.annotation.df) == "annotation"] <- "Best.Annotation.with.Analogues"
names(final.annotation.df)[names(final.annotation.df) == "smiles"] <- "Best.Annotation.with.Analogues.Smiles"

names(cytoscape.annotations) <- c("shared.name", "library_compound_name2")

cytoscape.annotations$library_compound_name2 <- ifelse(
  is.na(cytoscape.annotations$library_compound_name2) | cytoscape.annotations$library_compound_name2 == "NA", 
  "", 
  cytoscape.annotations$library_compound_name2
)
names(cytoscape.annotations) <- c("shared.name", "library.compound.name.2", "smiles", "confidence.level")

cytoscape <- cytoscape %>%
  full_join(cytoscape.annotations, by = "shared.name")
write.csv(cytoscape, "cytoscape-v2.csv")

samples.df <- final.annotation.df %>%
  select(feature.ID, feature.usi, Samples) %>%
  separate_rows(Samples, sep = "; ")
write.csv(samples.df, "samples-df.csv")

#More confidence.level assignments before finalising file.
for (i in 1:nrow(final.annotation.df)) {
  if (is.na(final.annotation.df$confidence.level[i]) && 
      (!is.na(final.annotation.df$ms2query.NPC.pathway[i]) || !is.na(final.annotation.df$canopus.NPC.pathway[i])) &&
      (!is.na(final.annotation.df$ms2query.NPC.pathway[i]) && final.annotation.df$ms2query.NPC.pathway[i] != "None")) { # Check if NOT NA before comparing
    final.annotation.df$confidence.level[i] <- "4"
  }
}

for (i in 1:nrow(final.annotation.df)) {
  if (is.na(final.annotation.df$confidence.level[i])) {
    final.annotation.df$confidence.level[i] <- "5"
  }
}

write.csv(final.annotation.df, "final-annotation-df.csv")