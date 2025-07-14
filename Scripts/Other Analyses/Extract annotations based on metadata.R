##Annotations across small-scale and large-scale

# Load necessary libraries
library(dplyr)
library(stringr)

#Dataset annotations
d125 <- read.csv("Y:/MA_BPA_Microbiome/Dataset-Annotations/HGMD_0125.csv") %>%   ##Tidying due to advanced script
  select(-gnps.MolecularFormula, -gnps.MonoisotopicMass, -gnps.ID.prob, -"gnps.compound.MolecularFormula", -"gnps.compound.MonoisotopicMass", 
         -"gnps.CID", -"authentic.standard.CID", -"authentic.standard.MolecularFormula",-"authentic.standard.MonoisotopicMass")
d108 <- read.csv("Y:/MA_BPA_Microbiome/Dataset-Annotations/HGMD_0108.csv")
d75 <- read.csv("Y:/MA_BPA_Microbiome/Dataset-Annotations/HGMD_0075.csv") %>%
  mutate(metadata = "Phe-Hex_LS")
d71 <- read.csv("Y:/MA_BPA_Microbiome/Dataset-Annotations/HGMD_0071.csv") %>%
  mutate(metadata = "Phe-Hex_LS")
d80 <- read.csv("Y:/MA_BPA_Microbiome/Dataset-Annotations/HGMD_0080.csv") %>%
  mutate(metadata = "Phe-Hex_LS")

#Metadata
metadata <- read.csv("Y:/MA_BPA_Microbiome/LCMS data/HGMD_0130_SS-Vs_LS[HILIC and Phe-Hex]/metadata.csv")
m108 <- filter(metadata, Dataset == "HGMD_0108") %>%
  filter(SS == 1)
m125 <- filter(metadata, Dataset == "HGMD_0125")
m125.ss <- filter(m125, SS == 1)
m125.ls <- filter(m125, LS == 1)

#Extract 108 and duplicate
search_pattern <- paste(unique(m108$filename), collapse = "|")
filtered_d108 <- d108 %>%
  filter(str_detect(Samples, search_pattern))
d108.ss <- filtered_d108 %>%
  mutate(metadata = "Phe-Hex_SS")
d108.ls <- filtered_d108 %>%
  mutate(metadata = "Phe-Hex_LS")

#Extract 125 for LS and SS separately
search_pattern <- paste(unique(m125.ss$filename), collapse = "|")
d125.ss <- d125 %>%
  filter(str_detect(Samples, search_pattern)) %>%
  mutate(metadata = "HILIC_SS")
search_pattern <- paste(unique(m125.ls$filename), collapse = "|")
d125.ls <- d125 %>%
  filter(str_detect(Samples, search_pattern)) %>%
  mutate(metadata = "HILIC_LS")

#Final combining
combined <- rbind(d75, d71) %>%
  rbind(d80) %>%
  rbind(d108.ss) %>%
  rbind(d108.ls) %>%
  rbind(d125.ss) %>%
  rbind(d125.ls)