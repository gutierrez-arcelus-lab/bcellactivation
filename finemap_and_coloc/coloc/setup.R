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

# Select genes with TSS within 250kb from the sentinel
bed38 <- 
    read_tsv(bed38_file, col_names = names(bed19)) |> 
    mutate(start = start + 2.5e5, end = end - 2.5e5)

gencode <- 
    "/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens/gencode.v39.primary_assembly.annotation.gtf.gz" |>
    read_tsv(comment = "#", col_names = FALSE, col_types = "ccciicccc") |>
    filter(X3 == "gene", X1 %in% paste0("chr", c(1:22, "X"))) |>
    mutate(gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
	   gene_id = str_remove(gene_id, "\\.\\d+$"),
	   gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+"),
	   tss = case_when(X7 == "+" ~ X4,
			   X7 == "-" ~ X5,
			   TRUE ~ NA)) |>
    select(chrom = 1, tss, gene_id, gene_name)

inner_join(gencode, bed38, join_by(chrom, between(tss, start, end))) |>
    select(locus, chrom, gene_id, gene_name) |>
    write_tsv(glue("./data/{gwas}_genes.tsv"))

# Set for coloc.susie
# Susie files from eQTL Catalogue only have position in GRCh38 and alleles, not rsid
# Bentham GWAS summary stats in GRCh38 to get positions

if (gwas == "Bentham") {
    
    gwas_38 <-
	"/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Bentham/harmonized/26502338-GCST003156-EFO_0002690.h.tsv.gz" |>
	read_tsv() |>
	select(rsid = variant_id, chrom = chromosome, pos = base_pair_location) |> 
	arrange(chrom, pos)
}

write_tsv(gwas_38, glue("./data/{gwas}_hg38_snps.tsv"))

# Slice 1Mb-windows into 100kb-windows to circumvent bug in eQTL Catalogue
# Filter by gene in the API returns wrong data
# And requesting data in 1Mb-windows sometimes go into the 100,000-row limit
windows_df <- 
    read_tsv(bed38_file, col_names = c("chrom", "start", "end", "locus")) |>
    mutate(chrom = str_remove(chrom, "^chr"))

windows_100kb <- 
    windows_df |>
    mutate(rg = map2(start, end, ~.x:.y)) |>
    select(locus, chrom, rg) |>
    unnest(cols = rg) |>
    group_by(locus) |>
    mutate(i = ntile(rg, 10)) |>
    group_by(locus, chrom, i) |>
    summarise(start = min(rg), end = max(rg)) |>
    ungroup() |>
    arrange(chrom, start, end)

write_tsv(windows_100kb, glue("./data/{gwas}_windows_100kb.tsv"))
