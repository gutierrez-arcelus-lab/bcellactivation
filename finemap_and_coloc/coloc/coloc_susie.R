library(tidyverse)

gwas_windows <-
    "../finemap/Bentham/data/windows.tsv" |>
    read_tsv(col_names = c("locus", "coords")) |>
    rowid_to_column()

gwas_hg38 <- 
    "../finemap/Bentham/data/summary_stats_hg38.tsv" |>
    read_tsv()

gwas_snp_df <- 
    gwas_hg38 |>
    mutate(variant = glue::glue("chr{chrom}_{pos}_{ref}_{alt}"),
	   snp_id = glue::glue("{rsid}-{ref}-{alt}")) |>
    select(locus, snp_id, variant)

gwas_genes <-
    read_tsv("./data/Bentham_genes.tsv")

gwas_lbf <- 
    glue::glue("../finemap/Bentham/data/susie/lbf/window{gwas_windows$rowid}.tsv") |>
    setNames(gwas_windows$locus) |>
    map_dfr(read_tsv, .id = "locus") |>
    inner_join(gwas_snp_df, join_by(locus, snp_id)) |>
    select(locus, variant, L1:L5)  

#    |>
#    {function(x) split(x, x$locus)}() |>
#    map(~select(., -locus)) |>
#    map(~pivot_longer(., L1:L5, names_to = "signal") |> 
#	pivot_wider(names_from = variant, values_from = value))

eqtl <- 
    "/lab-share/IM-Gutierrez-e2/Public/External_datasets/eQTLcatalogue/susie/QTS000001/QTD000002/QTD000002.lbf_variable.txt.gz" |>
    read_tsv()

eqtl |>
    filter(molecular_trait_id %in% gwas_genes$gene_id)


# subset for genes previously selected by TSS criterion
