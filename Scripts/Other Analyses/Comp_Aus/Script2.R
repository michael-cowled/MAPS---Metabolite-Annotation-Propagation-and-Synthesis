comp_list <- read.csv("comp_aust_list.csv")

# Assuming comp_list is your dataframe and Target is the column you want to check
comp_list$Metabolite_Present <- grepl("metabolite", comp_list$Target, ignore.case = TRUE)

# This will add a new column Metabolite_Present to comp_list indicating whether "metabolite" is present in each element of the Target column


# Assuming comp_list is your dataframe
colnames(comp_list)[colnames(comp_list) == "CAS.Number"] <- "cas"

# This will rename the column "Cas.Number" to "cas" in your dataframe comp_list


# Assuming cid_df and comp_list are your dataframes
merged_df <- merge(cid_df, comp_list, by = "cas", all = TRUE)

# This will merge the dataframes cid_df and comp_list by the common column "cas"
# The argument 'all = TRUE' ensures that all rows from both dataframes are retained in the merged dataframe, even if there are no matches in the other dataframe.


