library(rvest)
library(purrr)
library(stringr) # for str_trim

get_pubchem_titles <- function(cids) {
  
  get_title_single <- function(cid) {
    url <- paste0("https://pubchem.ncbi.nlm.nih.gov/compound/", cid)
    
    # Safely read the HTML page
    page <- tryCatch(
      rvest::read_html(url), 
      error = function(e) {
        warning(paste("Failed to retrieve page for CID:", cid, "- Error:", e$message))
        return(NA_character_)
      }
    )
    
    if (is.na(page) || !inherits(page, "xml_document")) {
      return(NA_character_)
    }
    
    # Extract the title tag content
    title_node <- rvest::html_node(page, "title")
    
    if (is.null(title_node)) {
      return(NA_character_)
    }
    
    title <- rvest::html_text(title_node, trim = TRUE)
    
    # Clean up the title string
    title <- gsub(" - PubChem", "", title, fixed = TRUE)
    title <- stringr::str_split(title, "\\|")[[1]][1]
    title <- stringr::str_trim(title)
    
    # Respect PubChem's rate limit
    Sys.sleep(0.2)
    
    return(title)
  }
  
  # Use purrr::map_chr to apply the function to each CID
  titles <- purrr::map_chr(cids, get_title_single)
  
  return(titles)
}



##Testing

my_data <- data.frame(CID = distinct_cids)

# Select the first three CIDs from the 'CID' column
first_three_cids <- my_data$CID[1:3]

# Now you can use 'first_three_cids' in your function call
# get_pubchem_titles(first_three_cids)

# Call the function on this smaller list
titles_for_test <- get_pubchem_titles(first_three_cids)

# Create a data frame to view the results
test_results <- data.frame(
  CID = first_three_cids,
  Title = titles_for_test
)

# Print the results
print(test_results)

##Actual

# 1. Get all unique CIDs from your main data frame
# Make sure to handle the case where 'my_data' is a data frame with a 'CID' column
distinct_cids <- unique(my_data$CID)

# 2. Call the function to get titles for all unique CIDs
# This will return a data frame with two columns: CID and Title
all_titles_df <- get_pubchem_titles(distinct_cids)

# 3. Join the new titles back to your original data frame
my_data_with_titles <- my_data %>%
  dplyr::left_join(all_titles_df, by = "CID")



# Now, 'my_data_with_titles' will contain a new column named 'Title'
# with the corresponding PubChem titles for each CID.

