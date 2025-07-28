##Query to find which features are in a specific fraction BUT not in the crude
#Use final.annotation.df2 so that best annotations are listed

##!!USER_DEFINED!!##
keep_sample <- "HGMA_1389" ##Sample of Interest
remove_sample <- "HGMA_1594" ##The crude to be removed
remove_sample2 <- "HGMA_1388" ##Other sample if needed

# Filter out rows containing the specified sample
final.annotation.df3 <- final.annotation.df2 %>%
  filter(str_detect(Samples, keep_sample)) %>%
  filter(!str_detect(Samples, remove_sample)) %>%
  filter(!str_detect(Samples, remove_sample2))





##TEST COMPARISON

keep_sample <- "HGMA_1388" ##Sample of Interest


# Filter out rows containing the specified sample
final.annotation.df4 <- final.annotation.df2 %>%
  filter(str_detect(Samples, keep_sample)) %>%
  filter(!str_detect(Samples, remove_sample))
