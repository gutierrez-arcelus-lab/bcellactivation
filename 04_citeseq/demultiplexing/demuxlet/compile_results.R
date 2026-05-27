library(tidyverse)

read_demuxlet <- function(f) {
    read_tsv(f) |>
    select(barcode = BARCODE, best = BEST) |>
    extract(best, c("status", "sample"), "([^-]+)-(.+)")
}

batches <- c("1984", "1988", "1990")

donors <- read_lines("../genotypes/data/mgb_sample_ids.txt")

demuxlet_df <- 
    glue::glue("./results/{batches}/results.best") |>
    setNames(batches) |>
    map_dfr(read_demuxlet, .id = "batch") |>
    select(barcode, batch, status, donor_id = sample) |>
    mutate(donor_id = case_when(status != "SNG" ~ status,
				TRUE ~ donor_id),
	   donor_id = case_when(status == "SNG" ~ str_extract(donor_id, "(\\d+)$"),
				TRUE ~ donor_id))

write_tsv(demuxlet_df, "./results/demuxlet_calls.tsv") 
