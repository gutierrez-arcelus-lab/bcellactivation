library(tidyverse)

gwas <- read_tsv("../data/gwas_catalog_v1.0.2.tsv", col_types = c(.default = "c"))

sle_gwas <- gwas %>%
    mutate(`P-VALUE` = as.numeric(`P-VALUE`)) %>%	
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
	   snp_current = SNP_ID_CURRENT,
	   p = `P-VALUE`) %>%
    distinct(chr, pos, .keep_all = TRUE)

sle_gwas %>%
    mutate(gene = ifelse(is.na(reported_gene), mapped_gene, reported_gene),
	   gene = na_if(gene, "NR")) %>%
    filter(!is.na(gene)) %>%
    select(-reported_gene, -mapped_gene) %>%
    separate_rows(gene, sep = ",") %>%
    separate_rows(gene, sep = " - ") %>%
    mutate(gene = trimws(gene),
	   gene = na_if(gene, "NA")) %>%
    filter(!is.na(gene)) %>%
    group_by(author, trait, gene) %>%
    slice(which.min(p)) %>%
    ungroup() %>% 
    group_by(gene) %>%
    summarise(author = paste(author, collapse = ","),
	      snp = paste(snp, collapse = ","),
	      p = paste(p, collapse = ",")) %>%
    ungroup() %>%
    write_tsv("./reported_genes.tsv")
