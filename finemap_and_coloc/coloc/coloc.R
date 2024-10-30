library(tidyverse)
library(glue)
library(coloc)

Sys.setenv(VROOM_CONNECTION_SIZE = 500000L)

run_coloc <- function(min_df) {

    eqtl_dataset <- list(beta = min_df$beta_qtl,
			 varbeta = min_df$se_qtl^2,
			 N = min_df$an_qtl/2,
			 MAF = min_df$maf_qtl,
			 type = "quant",
			 snp = min_df$variant)

    gwas_dataset <- list(beta = min_df$beta_gwas,
			 varbeta = min_df$se_gwas^2,
			 type = "cc",
			 snp = min_df$variant)

    coloc_res <- coloc.abf(dataset1 = eqtl_dataset, 
			   dataset2 = gwas_dataset)

    coloc_res
}

main <- function(qtl_file) {
    
    qtl_raw <- read_rds(qtl_file)

    qtl_data <- qtl_raw |>
	filter(map_lgl(result, ~!is.null(.))) |>
	select(-gene_id, -error) |>
	unnest(cols = result)

    if ( nrow(qtl_data) > 0 ) {	
    
	qtl_data <- qtl_data |>
	    group_by(molecular_trait_id) |>
	    filter(any(pvalue < 1e-5)) |>
	    ungroup()
    }

    if ( nrow(qtl_data) == 0 ) {
	
	out <- 
	    tibble(locus = NA, gene_name = NA, molecular_trait_id = NA, nsnps = NA, 
		   PP.H0.abf = NA, PP.H1.abf = NA, PP.H2.abf = NA, PP.H3.abf = NA, PP.H4.abf = NA)

	return(out)
    }

    min_df <- 
	inner_join(qtl_data, gwas_data, 
		   join_by(locus, variant), 
		   suffix = c("_qtl", "_gwas")) |>
	distinct(locus, gene_name, molecular_trait_id, 
		 variant, beta_qtl, beta_gwas, se_qtl, se_gwas, an_qtl = an, maf_qtl = maf)

    res_df <- 
	min_df |>
	group_by(locus, gene_name, molecular_trait_id) |>
	nest() |>
	ungroup() |>
	mutate(result = map(data, run_coloc))

    summary_df <- 
	res_df |>
	mutate(summ = map(result, "summary"),
	       summ = map(summ, ~enframe(.) |> pivot_wider(names_from = name, values_from = value))) |>
	select(locus, gene_name, molecular_trait_id, summ) |>
	unnest(cols = summ)

    summary_df
}

# GWAS data
gwas_38 <- 
    "./data/Bentham_hg38_snps.tsv" |>
    read_tsv()

gwas_data <- 
    "../finemap/Bentham/data/summary_stats.tsv" |>
    read_tsv() |>
    left_join(gwas_38, join_by(chrom, rsid)) |>
    mutate(variant = glue("chr{chrom}_{pos.y}_{ref}_{alt}")) |>
    select(locus, variant, beta, se, logp)

# QTL data
datasets <- read_tsv("./data/qtl_datasets.tsv")

qtl_files <- 
    glue("./data/qtls/Bentham_{datasets$dataset_id}.rds") |>
    setNames(datasets$dataset_id)

# Run
out_df <- map_dfr(qtl_files, main, .id = "dataset_id")
write_tsv(out_df, "./data/Bentham_coloc_results.tsv")
