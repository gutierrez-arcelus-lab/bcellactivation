library(tidyverse)

gwas <- read_tsv("../data/gwas_catalog_v1.0.2.tsv")

sle_gwas <- gwas %>%
    filter(`DISEASE/TRAIT` %in% c("Systemic lupus erythematosus",
				  "Lupus nephritis in systemic lupus erythematosus",
				  "Cutaneous lupus erythematosus",
				  "Childhood onset systemic lupus erythematosus")) %>%
    filter(`P-VALUE` < 5 * 1e-8) %>%
    select(author = `FIRST AUTHOR`,
	   trait = `DISEASE/TRAIT`,
	   region = REGION,
	   chr = CHR_ID,
	   pos = CHR_POS,
	   reported_gene = `REPORTED GENE(S)`,
	   mapped_gene = MAPPED_GENE,
	   snp = SNPS, 
	   snp_current = SNP_ID_CURRENT) %>%
    distinct(chr, pos, .keep_all = TRUE)

sle_gwas %>%
    separate_rows(reported_gene, sep = ",") %>%
    mutate(reported_gene = trimws(reported_gene)) %>%
    distinct(reported_gene) %>%
    pull(reported_gene) %>%
    write_lines("./reported_genes.tsv")
