library(tidyverse)
library(glue)

# Function to query eQTL Catalogue
request_associations <- function(dataset_id, chromosome_id, range_start, range_end) {
    
    page_size <- 1000
    page_start <- 0

    while (TRUE) {
	
	URL <- 
	    "https://www.ebi.ac.uk/eqtl/api/v2/datasets/{dataset_id}/associations?size={page_size}&start={page_start}&pos={chromosome_id}:{range_start}-{range_end}" |>
	    glue::glue()
	    
	r <- httr::GET(URL, httr::accept_json())

	# If Error == 500, try again up to 30 times
	ctn <- 0
	while (httr::status_code(r) == 500 && ctn < 30) {
	    r <- httr::GET(URL, httr::accept_json())
	    ctn <- ctn + 1L
	    Sys.sleep(5)
	}

	cont <- httr::content(r, "text", encoding = "UTF-8")

	# If the request was unsuccessful
	final_status <- httr::status_code(r) 
	if (final_status != 200) {
	    
	    #If we get no results at all, print error
	    if (page_start == 0) {
		stop(glue::glue("{final_status}\n{cont}"))
	    }
	    
	    #else just break
	    break
	}

	cont_df <- jsonlite::fromJSON(cont)
	
	if (page_start == 0) {
	    responses <- cont_df
	} else {
	    responses <- dplyr::bind_rows(responses, cont_df)
	}

	page_start <- page_start + page_size
    }
    
    return(responses)
}

safe_request <- safely(.f = request_associations)

# QTL dataset ID
dataset_id <- commandArgs(TRUE)[1]

# GWAS regions
gwas <- "Bentham"

windows_100kb <- read_tsv("./data/Bentham_windows_100kb.tsv") 

# Run query
res <-
    windows_100kb |>
    mutate(data = pmap(list(chrom, start, end),
		       function(x, y, z)
		       safe_request(dataset_id = dataset_id, 
				    chromosome_id = x,
				    range_start = y,
				    range_end = z)))

out <- 
    res |>
    mutate(result = map(data, "result"),
	   error = map(data, "error")) |>
    select(locus, result, error)

write_rds(out, glue("./data/qtls/{gwas}_{dataset_id}.rds"))
