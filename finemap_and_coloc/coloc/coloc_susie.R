library(tidyverse)
library(glue)
library(coloc)

# GRCh38 positions for GWAS variants
gwas_38 <- read_tsv("./data/Bentham_hg38_snps.tsv")

# GWAS windows
gwas_windows <-
    "../finemap/Bentham/data/windows.tsv" |>
    read_tsv(col_names = c("locus", "coords")) |>
    rowid_to_column() |>
    filter(locus == "IRF5")

# Genes with TSS within +- 250kb from lead variant
gwas_genes <-
    read_tsv("./data/Bentham_genes.tsv") |>
    filter(locus == "IRF5")

# Susie LBF data for GWAS
gwas_lbf <- 
    glue("../finemap/Bentham/data/susie/lbf/window{gwas_windows$rowid}.tsv") |>
    setNames(gwas_windows$locus) |>
    map_dfr(read_tsv, .id = "locus") |>
    separate(snp_id, c("rsid", "ref", "alt"), sep = "-") |>
    left_join(gwas_38, join_by(rsid)) |>
    mutate(variant = glue("chr{chrom}_{pos}_{ref}_{alt}")) |>
    select(locus, variant, rsid, L1:L5)

# Susie LBF data for QTL (header)
qtl_names <- 
    "/lab-share/IM-Gutierrez-e2/Public/External_datasets/eQTLcatalogue/susie/QTS000001/QTD000001/QTD000001.lbf_variable.txt.gz" |>
    read_tsv(n_max = 1) |>
    colnames()

# eQTL Catalogue datasets
qtl_datasets <-
    read_tsv("./data/qtl_datasets.tsv") |>
    select(dataset_id, study_label, sample_group, tissue_label, condition_label, quant_method)


run_coloc <- function(dataset_id) {

    qtl <- 
	glue("./data/qtls_susie/Bentham/{dataset_id}.lbf_variable.txt.gz") |>
	read_tsv(col_names = qtl_names) |>
	mutate(gene_id = str_extract(molecular_trait_id, "ENSG\\d+")) |>
	select(gene_id, everything()) |>
	filter(gene_id %in% gwas_genes$gene_id)

    if (nrow(qtl) == 0) return("No QTL data")

    # Format LBF data
    gwas_lbf_matrix <-
	gwas_lbf |>
	filter(variant %in% qtl$variant) |>
	pivot_longer(L1:L5, names_to = "cs") |>
	select(cs, variant, value) |>
	pivot_wider(names_from = variant, values_from = value) |>
	column_to_rownames("cs") |>
	as.matrix()

    qtl_lbf_matrix <- 
	qtl |>
	{function(x) split(x, x$molecular_trait_id)}() |>
	map(~pivot_longer(., lbf_variable1:lbf_variable10, names_to = "cs") |>
	    select(cs, variant, value) |>
	    filter(variant %in% colnames(gwas_lbf_matrix)) |>
	    pivot_wider(names_from = variant, values_from = value) |>
	    column_to_rownames("cs") |>
	    select(all_of(colnames(gwas_lbf_matrix))) |>
	    as.matrix())

    # Run coloc
    map_dfr(qtl_lbf_matrix, 
	    ~coloc.bf_bf(gwas_lbf_matrix, .)$summary, 
	    .id = "molecular_trait_id") |>
    as_tibble()
}

res_df <-
    qtl_datasets |>
    mutate(results = map(dataset_id, run_coloc))

out_df <- 
    res_df |>
    filter(map_lgl(res_df$results, is.data.frame)) |>
    unnest(cols = results)

write_tsv(out_df, "./data/coloc_susie/results.tsv")


# Filter results

## GWAS
gwas_stats <-
    "../finemap/Bentham/data/summary_stats.tsv" |>
    read_tsv() |>
    select(-chrom, -pos) |>
    left_join(gwas_38, join_by(rsid)) |>
    mutate(variant = glue("chr{chrom}_{pos}_{ref}_{alt}")) |>
    distinct(variant, logp)

out_gwas_filtered <- 
    out_df |>
    left_join(gwas_stats, join_by(hit1 == variant)) |>
    filter(logp > 5) |>
    select(-logp)

## QTL
hits <- 
    out_gwas_filtered |>
    distinct(dataset_id, molecular_trait_id, hit2)

qtl_files <- 
    glue("./data/qtls/Bentham_{qtl_datasets$dataset_id}.rds") |>
    setNames(qtl_datasets$dataset_id)

qtl_data <- 
    qtl_files |>
    map_dfr(function(x) {
		dat <- 
		    read_rds(x) |>
		    filter(map_lgl(result, ~!is.null(.))) |>
		    select(-gene_id, -error) |>
		    unnest(cols = result)

		if ( nrow(dat) == 0 ) {

		    out <- tibble(molecular_trait_id = NA, variant = NA, nlog10p = NA)

		    return(out)
		}
	    
		dat |>   
		filter(variant %in% hits$hit2) |>
		select(molecular_trait_id, variant, nlog10p)
	    },
	    .id = "dataset_id")

out_filtered <- 
    out_gwas_filtered |>
    left_join(qtl_data, join_by(dataset_id, molecular_trait_id, hit2 == variant)) |>
    filter(nlog10p > 5) |>
    select(-nlog10p)
   
out_filtered |>
    mutate(gene_id = str_extract(molecular_trait_id, "ENSG\\d+")) |>
    inner_join(gwas_genes, join_by(gene_id)) |>
    select(dataset_id:quant_method, gene_id, gene_name, molecular_trait_id, nsnps, idx1, idx2, hit1, hit2, 
	   PP.H0.abf:PP.H4.abf) |>
    arrange(desc(PP.H4.abf)) |>	
    write_tsv("./data/coloc_susie/results_filtered.tsv")
