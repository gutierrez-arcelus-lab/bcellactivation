library(coloc)
library(GenomicRanges)
library(seqminer)
library(tidyverse)
library(furrr)

# functions
fetch <- function(f, r) {
    tabix.read.table(tabixFile = f, 
		     tabixRange = r) |>
    as_tibble()
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
    
    eqtl <- fetch(f = ftp, r = region) |>
	as_tibble() |>
	setNames(header) |>
	filter(gene_id %in% genes)

    eqtl_filtered <- eqtl |>
	mutate(rsid = sub("\\\r$", "", rsid)) |>
	filter(rsid != "NA") |>
	group_by(molecular_trait_id, gene_id) |>
	filter(any(pvalue < 5e-5)) |>
	ungroup() |>
	add_count(molecular_trait_id, gene_id, rsid) |>
	filter(n == 1) |>
	select(-n) |>
	add_count(molecular_trait_id, gene_id, type, position) |>
	filter(n == 1) |>
	select(-n)

    if ( nrow(eqtl_filtered) == 0 ) {
	return(tibble(gene_id = NA, molecular_trait_id = NA, nsnps = NA, 
		      h0 = NA, h1 = NA, h2 = NA, h3 = NA, h4 = NA))
    }

    eqtl_cat_df <- eqtl_filtered |>
	mutate(eqtl_sample_size = an/2L,
	       eqtl_varbeta = se^2) |>
	select(gene_id, molecular_trait_id, rsid, position, 
	       eqtl_sample_size, eqtl_maf = maf, eqtl_beta = beta, eqtl_varbeta)

    min_df <- inner_join(eqtl_cat_df, gwas_sum, by = c("rsid", "position" = "pos")) |>
	filter(!is.na(eqtl_varbeta)) |>
	filter(gwas_beta != 0, gwas_varbeta != 0)

    coloc_results <- min_df |>
	unite("id", c(gene_id, molecular_trait_id), sep = ",") |>
	(\(x) split(x, x$id))() |>
	map(run_coloc)

    coloc_results_summary <- map(coloc_results, "summary") |>
	map_dfr(function(x) as_tibble(x, rownames = "stat"), .id = "id") |>
	pivot_wider(names_from = stat, values_from = value) |>
	select(id, nsnps, h0 = PP.H0.abf, h1 = PP.H1.abf, h2 = PP.H2.abf, h3 = PP.H3.abf, h4 = PP.H4.abf) |>
	separate(id, c("gene_id", "molecular_trait_id"), sep = ",") |>
	select(gene_id, molecular_trait_id, everything())

    coloc_results_summary
}


# Cmd arguments
pargs <- commandArgs(TRUE)[1:2]
dataset <- pargs[1]
region_id <- pargs[2]

# eQTL Catalogue column names
column_names <- read_lines("eqtl_column_names.txt")

# eQTL catalogue URLs
ftps <- "data/coloc_input/ftp_urls.txt" |>
    read_lines()

study_ids <- basename(ftps) |> 
    sub("\\.all\\.tsv\\.gz", "", x = _)

names(ftps) <- study_ids

# GWAS summ stats
gwas_data <- "data/coloc_input/gwas_data_%s.tsv" |>
    sprintf(dataset) |>
    read_tsv() |>
    group_nest(region, .key = "gwas_sum")

# regions
region_df <- "data/coloc_input/regions_%s.tsv" |>
    sprintf(dataset) |>
    read_tsv()

region_df_i <- region_df |>
    filter(region == region_id) |>
    left_join(gwas_data, by = "region")

# genes to consider (TSS 250kb from GWAS hit)
genes_dat <- "data/coloc_input/gene_annots_%s.tsv" |>
    sprintf(dataset) |>
    read_tsv() |>
    filter(region == region_id) |>
    select(gene_id, gene_name)

# Set up parallel architecture
plan(cluster, workers = length(availableWorkers()))

# Run analysis
analysis_df <- enframe(ftps) |>
    setNames(c("study", "ftp"))

coloc_res_list <- list()
error_i <- seq_len(nrow(analysis_df))
i <- 1

while ( length(error_i) > 0 && i <= 30) {

    coloc_res_list[[i]] <- analysis_df |>
	slice(error_i) |>
	mutate(res = future_map(ftp,  
				safely(function(x) main(ftp = x, 
							region = region_df_i$coord, 
							gwas_sum = region_df_i$gwas_sum[[1]],
							genes = genes_dat$gene_id,
							header = column_names))))

    errors <- coloc_res_list[[i]] |>
	mutate(err = map(res, "error")) |>
	pull(err)

    error_i <- errors |>
	map(~!is.null(.)) |>
	unlist() |>
	which()

    i <- i + 1
}

if ( length(error_i) > 0) stop("Could not retrieve or process eQTL catalogue data.")

results <- coloc_res_list |>
    bind_rows() |>
    mutate(dat = map(res, "result")) |>
    select(study, dat) |>
    unnest(cols = dat) |>
    left_join(genes_dat, by = "gene_id") |>
    select(study, gene_id, gene_name, everything())

out_name <- file.path("results", 
		      sprintf("%s_region%s.tsv", dataset, region_id))

write_tsv(results, out_name)
