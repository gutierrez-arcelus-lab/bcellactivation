library(tidyverse)
library(httr)
library(rvest)
library(jsonlite)
library(glue)


if ( !file.exists("data") ) dir.create("data")

# eQTL Catalogue
URL <- "https://www.ebi.ac.uk/eqtl/api/v2/datasets/?size=1000"

## Make a request
r <- GET(URL, accept_json())

## Check status
status_code(r)

## Extract content and convert content to dataframe
datasets <- 
    content(r, "text", encoding = "UTF-8") |>
    fromJSON() |>
    as_tibble()

write_tsv(datasets, "./data/qtl_datasets.tsv")


# liftOver GWAS windows
gwas <- "Bentham"
bed19_file <- glue("./data/{gwas}_hg19.bed")
bed38_file <- glue("./data/{gwas}_hg38.bed")
fail_file <- glue("./data/{gwas}.failTolift.txt")
chain_file <- "/reference_databases/ReferenceGenome/liftover_chain/hg19/hg19ToHg38.over.chain.gz"

bed19 <-
    glue("../finemap/{gwas}/data/windows.tsv") |>
    read_tsv(col_names = c("locus", "coords")) |>
    extract(coords, c("chrom", "start", "end"), "(.+):(.+)-(.+)") |>
    mutate(chrom = paste0("chr", chrom)) |>
    select(chrom, start, end, locus)

write_tsv(bed19, bed19_file, col_names = FALSE)

command <- sprintf("liftOver %s %s %s %s", bed19_file, chain_file, bed38_file, fail_file)
system(command)
