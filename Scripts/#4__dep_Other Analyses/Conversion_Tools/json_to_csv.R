library(jsonlite)
json_data <- fromJSON("metabolite.index.extra.json.gz")
json_df <- as.data.frame(json_data)

json_df$HMDB_ID <- as.character(json_df$HMDB_ID)
json_df$Chemical_Name <- as.character(json_df$Chemical_Name)

write.csv(json_df, "HMDB.csv")