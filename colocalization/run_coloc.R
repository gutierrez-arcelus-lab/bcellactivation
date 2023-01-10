library(coloc)
library(GenomicRanges)
library(seqminer)
library(tidyverse)
library(furrr)

# functions
import_eqtl <- function(ftp_path, region, header, genes) {
 
    fetch <- safely(function(f, r) {
	tabix.read.table(tabixFile = f, 
			 tabixRange = r) |>
	as_tibble()
    })

    fetch_table <- fetch(f = ftp_path, r = region) 
    error <- fetch_table$error
	
    i <- 1
    while ( !is.null(error) & i < 30 ) {
	Sys.sleep(10)
	fetch_table <- fetch(f = ftp_path, r = region) 
	error <- fetch_table$error
	cat("Retrying...\n")
	i <- i + 1
    }
    
    if ( !is.null(fetch_table$error) ) stop("Can't retrieve dataset")

    out <- as_tibble(fetch_table$result)
    colnames(out) <- header

    filter(out, gene_id %in% genes)
}

run_coloc <- function(min_df) {

    eqtl_dataset <- list(beta = min_df$eqtl_beta,
			 varbeta = min_df$eqtl_varbeta,
			 N = min_df$eqtl_sample_size,
			 MAF = min_df$eqtl_maf,
			 type = "quant",
			 snp = min_df$rsid)

    gwas_dataset <- list(beta = min_df$gwas_beta,
			 varbeta = min_df$gwas_varbeta,
			 type = "cc",
			 snp = min_df$rsid)

    coloc_res <- coloc.abf(dataset1 = eqtl_dataset, 
			   dataset2 = gwas_dataset)

    coloc_res
}

main <- function(ftp, region, gwas_sum, genes, header) {
    
    eqtl <- import_eqtl(ftp, region, header, genes)

    eqtl_filtered <- eqtl %>%
	mutate(rsid = sub("\\\r$", "", rsid)) %>%
	filter(rsid != "NA") %>% 
	group_by(molecular_trait_id, gene_id) %>%
	filter(any(pvalue < 5e-5)) %>%
	ungroup() %>%
	add_count(molecular_trait_id, gene_id, rsid) %>%
	filter(n == 1) %>%
	select(-n) %>%
	add_count(molecular_trait_id, gene_id, type, position) %>%
	filter(n == 1) %>%
	select(-n)

    eqtl_cat_df <- eqtl_filtered %>%
	mutate(eqtl_sample_size = an/2L,
	       eqtl_varbeta = se^2) %>%
	select(gene_id, molecular_trait_id, rsid, position, 
	       eqtl_sample_size, eqtl_maf = maf, eqtl_beta = beta, eqtl_varbeta)

    min_df <- inner_join(eqtl_cat_df, gwas_sum, by = c("rsid", "position" = "pos")) %>%
	filter(!is.na(eqtl_varbeta)) %>%
	filter(gwas_beta != 0, gwas_varbeta != 0)

    coloc_results <- min_df %>%
	unite("id", c(gene_id, molecular_trait_id), sep = ",") %>%
	split(.$id) %>% 
	map(run_coloc)

    if (length(coloc_results) == 0) return(tibble())

    coloc_results_summary <- map(coloc_results, "summary") |>
	map_dfr(function(x) as_tibble(x, rownames = "stat"), .id = "id") |>
	pivot_wider(names_from = stat, values_from = value) |>
	select(id, nsnps, h0 = PP.H0.abf, h1 = PP.H1.abf, h2 = PP.H2.abf, h3 = PP.H3.abf, h4 = PP.H4.abf) |>
	separate(id, c("gene_id", "molecular_trait_id"), sep = ",") |>
	left_join(genes_df) |>
	select(gene_id, gene_name, molecular_trait_id, everything())

    coloc_results_summary
}

# Directoy
rundir <- "/lab-share/IM-Gutierrez-e2/Public/vitor/lupus/colocalization"

# Cmd arguments
pargs <- commandArgs(TRUE)[1:2]
array_id <- pargs[1]
procs <- pargs[2]

# eQTL Catalogue column names
column_names <- read_lines(file.path(rundir, "eqtl_column_names.txt"))

# GWAS summ stats
gwas_data <- file.path(rundir, "data/coloc_input/gwas_data.tsv") |>
    read_tsv() |>
    group_nest(region, .key = "gwas_sum")

# eQTL catalogue URLs
ftps <- file.path(rundir, "data/coloc_input/ftp_urls.txt") |>
    read_lines()

study_ids <- basename(ftps) |> 
    sub("\\.all\\.tsv\\.gz", "", x = _)

names(ftps) <- study_ids

# regions
region_df <- file.path(rundir, "data/coloc_input/region.tsv") |>
    read_tsv() |>
    filter(region == array_id) |>
    left_join(gwas_data, by = "region")

# genes to consider (TSS 250kb from GWAS hit)
genes_df <- file.path(rundir, "data/coloc_input/gene_annots.tsv") |>
    read_tsv()

# Run analysis
plan(multisession, workers = procs)

coloc_df <- 
    future_map_dfr(ftps, 
		   function(x) main(ftp = x, 
				    region = region_df$coord, 
				    gwas_sum = region_df$gwas_sum[[1]],
				    genes = genes_df$gene_id,
				    header = column_names),
		   .id = "study")

out_name <- file.path(rundir, "results", paste0("bentham_region", array_id, ".tsv"))

write_tsv(coloc_df, out_name)
