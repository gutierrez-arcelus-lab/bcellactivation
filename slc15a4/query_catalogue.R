library(tidyverse)
library(httr)
library(glue)
library(jsonlite)

# Function to request data in a region
request_associations <- function(dataset_id, chromosome_id, position, distance = 5e5) {
    
    page_size <- 1000
    page_start <- 0
    range_start <- position - distance
    range_end <- position + distance

    while (TRUE) {
	
	URL <- glue("https://www.ebi.ac.uk/eqtl/api/v2/datasets/{dataset_id}/associations?size={page_size}&start={page_start}&pos={chromosome_id}:{range_start}-{range_end}")
	r <- GET(URL, accept_json())

	# If Error == 500, repeat
	ctn <- 0
	while (status_code(r) == 500 && ctn < 30) {
	    r <- GET(URL, accept_json())
	    ctn <- ctn + 1L
	    Sys.sleep(5)
	    
	    #If we get no results at all, print error
	    if (page_start == 0 & ctn >= 30) {
		print(glue("Error {status_code(r)}"))
		print(cont)
		return()
	    }
	}

	cont <- content(r, "text", encoding = "UTF-8")

	# If the request was unsuccessful
	if (status_code(r) != 200) {
	    
	    #If we get no results at all, print error
	    if (page_start == 0) {
		print(glue("Error {status_code(r)}"))
		print(cont)
		return()
	    }
	    
	    #else just break
	    break
	}

	cont_df <- fromJSON(cont)
	
	if (page_start == 0) {
	    responses <- cont_df
	} else {
	    responses <- bind_rows(responses, cont_df)
	}

	page_start <- page_start + page_size
    }
    
    return(responses)
}

# Region of interest
chr <- 12
pos <- 128794319

# eQTL catalogue dataset
query_id <- commandArgs(TRUE)[1]

# Request eQTL catalogue data
associations <- 
    request_associations(dataset_id = query_id,
			 chromosome_id = chr,
			 position = pos) |>
    as_tibble()

if (is.data.frame(associations)) {
    
    out_file <- system("echo $TEMP_WORK", intern = TRUE) |>
	file.path("eqtl_catalogue/%s_chr%s_%d-%d.tsv") |>
	sprintf(query_id, chr, pos - 5e5L, pos + 5e5L)

    write_tsv(associations, out_file)
}

