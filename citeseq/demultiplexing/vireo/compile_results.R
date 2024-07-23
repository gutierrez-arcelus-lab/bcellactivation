library(tidyverse)

batches <- c("1984", "1988", "1990")

vireo_df <-
    glue::glue("./results/{batches}/donor_ids.tsv") |>
    setNames(batches) |>
    map_dfr(read_tsv, .id = "batch") |>
    mutate(status = case_when(donor_id %in% c("doublet", "unassigned") ~ donor_id, 
			      TRUE ~ "singlet"),
	   donor_id = case_when(status == "singlet" ~ str_extract(donor_id, "(\\d+)$"),
				TRUE ~ donor_id)) |>
    select(barcode = cell, batch, status, donor_id)

write_tsv(vireo_df, "./results/vireo_calls.tsv")
