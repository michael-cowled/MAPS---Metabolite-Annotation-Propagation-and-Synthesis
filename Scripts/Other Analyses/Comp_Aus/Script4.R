# Load the required library
library(VennDiagram)

setA <- merged_df2[merged_df2$Metabolite_Present >= 1, ]
setB <- merged_df2[merged_df2$presence.y >= 1, ]
setC <- merged_df2[merged_df2$Brunda.s.STD.mix.1.5 >= 1 & !is.na(merged_df2$Brunda.s.STD.mix.1.5), ]

# Create Venn diagram with increased font size and style for labels
venn.plot <- venn.diagram(
  x = list(A = rownames(setA), B = rownames(setB), C = rownames(setC)),
  category.names = c("Metabolite", "Natural Product", "Already an LCMS standard"),
  filename = NULL,
  cex = 2,  # Adjust font size (increase or decrease as needed)
  fontfamily = "sans",  # Set font family
  fontface = "bold",     # Set font style
  category.fontsize = 2
)

# Plot the Venn diagram
grid.draw(venn.plot)


##Which standards are in A and B but not C:
# Assuming setA, setB, and setC are defined as in your code

# Merge setA and setB by 'cas'
common_AB <- merge(setA, setB, by = "cas", all = TRUE, suffixes = c("_A", "_B"))

# Merge the result with setC by 'cas'
merged_all <- merge(common_AB, setC, by = "cas", all = TRUE)

# Get entries present only in setA and setB but not in setC
only_AB_not_C <- merged_all[is.na(merged_all$Metabolite_Present) &
                              is.na(merged_all$presence.y) &
                              !is.na(merged_all$Brunda.s.STD.mix.1.5), ]

# Display the resulting dataframe
print(only_AB_not_C)
