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

genes_df <- 
    glue("./data/{gwas}_genes.tsv") |>
    read_tsv()

run_df <-
    left_join(genes_df, windows_df, join_by(locus))

# Run query
res <-
    run_df |>
    mutate(data = pmap(.l = list(chrom, gene_id, start, end),
		       .f = function(chr, g, s, e)
			   safe_request(dataset_id = dataset_id, 
					chromosome_id = chr,
					gene_id = g,
					range_start = s,
					range_end = e)))

out <- 
    res |>
    mutate(result = map(data, "result"),
	   error = map(data, "error")) |>
    select(locus, gene_id, gene_name, result, error)

write_rds(out, glue("./data/qtls/{gwas}_{dataset_id}.rds"))
