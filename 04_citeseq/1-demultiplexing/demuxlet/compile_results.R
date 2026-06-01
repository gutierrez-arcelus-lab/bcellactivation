# ==============================================================================
# Description:  Parses and aggregates the raw demuxlet assignments (.best files) 
#               across all Cell Ranger library batches. It extracts the singlet/
#               doublet status and cleans the donor IDs to create a single master 
#               metadata table for downstream Seurat/single-cell filtering.
# Input:        1. ./results/{batch}/results.best (Demuxlet raw outputs)
# Output:       ./results/demuxlet_calls.tsv (Compiled cell assignments)
# ==============================================================================

library(tidyverse)
library(glue)

# ------------------------------------------------------------------------------
# 1. Function Definitions
# ------------------------------------------------------------------------------
# Custom function to read and parse the demuxlet output files

read_demuxlet <- function(f) {
    read_tsv(f) |>
    select(barcode = BARCODE, best = BEST) |>
    extract(best, c("status", "sample"), "([^-]+)-(.+)")
}

# ------------------------------------------------------------------------------
# 2. Data Parsing & Aggregation
# ------------------------------------------------------------------------------

# Define the library batches run through Cell Ranger
batches <- c("1984", "1988", "1990")

# Dynamically generate paths, read the files, and bind them into one dataframe
demuxlet_df <- 
    glue("./results/{batches}/results.best") |>
    setNames(batches) |>
    map_dfr(read_demuxlet, .id = "batch") |>
    select(barcode, batch, status, donor_id = sample) |>
    # Clean up the donor_id column for downstream processing:
    mutate(
	   # If the cell is a Doublet (DBL) or Ambiguous (AMB), overwrite the 
	   # complex multi-donor ID string with just the status (e.g., "DBL") 
	   # so they can be easily filtered out in Seurat.
	   donor_id = case_when(status != "SNG" ~ status,
				TRUE ~ donor_id),
	   # If the cell is a valid Singlet (SNG), extract just the trailing numeric 
	   # donor ID (removing any leading prefix text) for clean metadata.
	   donor_id = case_when(status == "SNG" ~ str_extract(donor_id, "(\\d+)$"),
				TRUE ~ donor_id))

# ------------------------------------------------------------------------------
# 4. Export Table
# ------------------------------------------------------------------------------
write_tsv(demuxlet_df, "./results/demuxlet_calls.tsv") 
