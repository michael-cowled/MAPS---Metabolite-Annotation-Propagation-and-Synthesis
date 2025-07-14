check_presence <- function(cid_df, data_url) {
  if (!require(httr)) {
    install.packages("httr")
    library(httr)
  }
  if (!require(jsonlite)) {
    install.packages("jsonlite")
    library(jsonlite)
  }
  
  presence_counts <- numeric(length(cid_df$cid))
  
  for (i in seq_along(cid_df$cid)) {
    url <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/", cid_df$cid[i], "/JSON")
    print(i)
    # Try fetching data, handle errors gracefully
    tryCatch({
      response <- GET(url)
      if (response$status_code == 200) {
        raw.data <- fromJSON(content(response, "text", encoding = "UTF-8"))
        presence_absence <- grepl("natural product", raw.data$Record$Section$Section, ignore.case = TRUE)
        presence_counts[i] <- sum(presence_absence)
      } else {
        presence_counts[i] <- NA
        warning(paste("Failed to fetch data for CID", cid_df$cid[i], "with status code:", response$status_code))
      }
    }, error = function(e) {
      presence_counts[i] <- NA
      warning(paste("Failed to fetch data for CID", cid_df$cid[i], "with error:", conditionMessage(e)))
    })
  }
  
  return(presence_counts)
}

# Usage example:
if (!file.exists("cas-to-cid.csv")) {
  stop("cas-to-cid.csv not found. Please provide the file.")
}

cid_df <- read.csv("cas-to-cid.csv")
presence_counts <- check_presence(cid_df, raw.data)
cid_df$presence <- presence_counts
print(presence_counts)