library(DBI)
cid_db_con <- dbConnect(RSQLite::SQLite(), "~/PubChem_Indexed.sqlite")

# Example CID list (replace with your actual vector of CIDs)
cid_list <- c(25163935)

# Convert the list to a string for SQL
cid_str <- paste(cid_list, collapse = ", ")

# Run SQL query – one row per CID
query <- sprintf("
  SELECT CID, Title, SMILES, IUPAC, Formula, \"Monoisotopic.Mass\"
  FROM pubchem_data
  WHERE CID IN (%s)
  GROUP BY CID
", cid_str)

# Execute query
results <- dbGetQuery(cid_db_con, query)

# View the results
print(results)


#Note: for instances where Title is NA, the IUPAC is used

# Close connection at end
if (dbIsValid(cid_db_con)) {
  dbDisconnect(cid_db_con)
  message("[DB DISCONNECT] Closed global DB connection.")
}