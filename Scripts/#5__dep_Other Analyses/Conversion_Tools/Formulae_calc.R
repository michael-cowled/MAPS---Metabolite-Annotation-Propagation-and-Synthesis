##Converts formulae to exact mass

library(MassTools)
library(dplyr)

formulae <- read.csv("Comp_Aust_List_with_SMILES_and_CID.csv")

i <- 8001
while (i <= 8437) {
  print(i)
  print(formulae[i,2])
  print(as.numeric(getExactMass(formulae[i, 10])))
  formulae[i,6] <- as.numeric(getExactMass(formulae[i, 10]))
  i <- i+1
}

cas <- select(formulae, -Column1)

###NOTE: will crash if formulae contain a "." like ".H2O", or is N/A

write.csv(formulae, "CompAust-extra-metadata-with-exact-mass_part9.csv")

#---#
##Converts Cas# to SMILES in prep for CFM-ID

#install.packages("remotes")
#remotes::install_github("harveyl888/chem")
library(chem)
library(tcltk)

cas<-cbind(cas,cas)
u <- 1:nrow(cas)
print(u)
pb <- tkProgressBar("progress","done %", 0, 100)
for (i in u) {
  cas[i,2]<-cas_to_smiles(cas[i,1])
  info <- sprintf("done %d%%", round(i*100/length(u)))
  setTkProgressBar(pb, i*100/length(u), sprintf("progress (%s)", info), info)
  print(i)
}

colnames(cas)[colnames(cas) == "X"] <- "SMILES"

write.csv(cas, "CompAust-extra-metadata-with-exact-mass-and-smiles.csv")



#### FOR complete Comp Aust

cas <- read.csv("Comp_Aust_CAS_only.csv")

write.csv(cas, "CompAust_to_smiles.csv")




##Add SMILE to whole Comp_Aust list

Comp_Aust_List <- read.csv("comp_aust_list.csv")
colnames(cas) <- c("Cas", "SMILES", "cas2", "Catalog.Number")
cas <- select(cas, -cas2)
# Assuming Comp_Aust_List and cas are your dataframes

# Merge the dataframes by Catalog.Number
merged_df <- merge(Comp_Aust_List, cas, by = "Catalog.Number")

write.csv(merged_df, "Comp_Aust_List_with_SMILES_and_CID.csv")
