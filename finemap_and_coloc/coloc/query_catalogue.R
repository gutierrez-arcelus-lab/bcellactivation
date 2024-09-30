library(tidyverse)
library(glue)

# Function to query eQTL Catalogue
source("./request_associations.R")
safe_request <- safely(.f = request_associations)

# QTL dataset ID
dataset_id <- commandArgs(TRUE)[1]

# GWAS regions
gwas <- "Bentham"

windows_df <- 
    glue("./data/{gwas}_hg38.bed") |>
    read_tsv(col_names = c("chrom", "start", "end", "locus")) |>
    mutate(chrom = str_remove(chrom, "^chr"))

# Run query
res <-
    windows_df |>
    mutate(data = pmap(.l = list(chrom, start, end),
		       .f = function(x, y, z)
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

