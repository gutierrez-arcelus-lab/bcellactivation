# Function to request data in a region
request_associations <- function(dataset_id, chromosome_id, gene_id, range_start, range_end) {
    
    page_size <- 1000
    page_start <- 0

    while (TRUE) {
	
	URL <- 
	    "https://www.ebi.ac.uk/eqtl/api/v2/datasets/{dataset_id}/associations?size={page_size}&start={page_start}&pos={chromosome_id}:{range_start}-{range_end}&gene_id={gene_id}" |>
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
	    responses <- bind_rows(responses, cont_df)
	}

	page_start <- page_start + page_size
    }
    
    return(responses)
}
