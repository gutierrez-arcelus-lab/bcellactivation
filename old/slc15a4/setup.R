library(tidyverse)
library(rvest)
library(httr)
library(glue)
library(jsonlite)

#dir.create("data")

bentham_window <- 
    tibble(chr = "12", pos = 129278864) |> 
    mutate(start = pos - 2.5e5, 
	   stop = pos + 2.5e5) |>
    select(chr, start, stop)

# Summary stats
bentham_stats <- 
    "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Bentham/opengwas/ebi-a-GCST003156.vcf.gz" |>
    read_tsv(comment = "##") |>
    janitor::clean_names() |>
    select(chr = number_chrom, pos, id, ref, alt, format, stats = 10) |>
    inner_join(bentham_window, join_by(chr, between(pos, start, stop))) |>
    separate_rows(format:stats, sep = ":") |>
    select(chr, pos, rsid = id, ref, alt, format, stats) |>
    pivot_wider(names_from = format, values_from = stats) |>
    select(chr, pos, rsid, ref, alt, beta = ES, se = SE, logp = LP) |>
    mutate_at(vars(beta, se, logp), as.numeric)

# Save summ stats
write_tsv(bentham_stats, "./data/summstats_bentham.tsv")

# Save variant rs IDs
bentham_stats |>
    pull(rsid) |>
    write_lines("./data/rsids_bentham.txt")

#get GRCh38 positions
system("./get_dbsnp.sh")

hg38 <- 
    "./data/gwas_rsid_dbsnp.vcf" |>
    read_tsv(comment = "##") |>
    janitor::clean_names()

hg38_window <- hg38 |>
    group_by(chr = number_chrom) |>
    summarise(start = min(pos), stop = max(pos)) |>
    ungroup()

write_tsv(hg38_window, "./data/window_bentham.tsv")

# eQTL Catalogue
#All datasets will be pulled if this parameter is bigger than the actual number of datasets
max_pulled_rows <- 1000 

URL <- glue("https://www.ebi.ac.uk/eqtl/api/v2/datasets/?size={max_pulled_rows}")

# Make a request
r <- GET(URL, accept_json())

# Check status
status_code(r)

# Extract content and convert content to dataframe
datasets <- 
    content(r, "text", encoding = "UTF-8") |>
    fromJSON() |>
    as_tibble()


write_tsv(datasets, "./data/eqtl_datasets.tsv")
