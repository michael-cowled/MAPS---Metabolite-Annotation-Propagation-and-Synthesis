dataset <- read.csv("Y:/MA_BPA_Microbiome/Dataset-Annotations/HGMD_0128.csv")

usi_counts_dataset <- dataset %>%
  group_by(confidence.level) %>%
  summarize(total_usi = n_distinct(feature.usi), .groups = "drop")

print(usi_counts_dataset)

