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

coords <- sle_gwas %>%
    select(chr, pos) %>%
    drop_na() %>%
    mutate(chr = factor(chr, levels = c(1:22, "X")),
	   pos = as.numeric(pos)) %>%
    arrange(chr, pos) %>%
    mutate(start = pos - 2e5 + 1L,
	   end = pos + 2e5 - 1L)

gencode <- read_tsv("../data/gencodev38_genes.bed", col_names = FALSE) %>%
    extract(X4, c("gene_id", "gene_name"), "([^_]+)_(.+)") %>%
    select(chr = X1, tss = X2, gene_id, gene_name) %>%
    mutate(chr = sub("chr", "", chr),
	   chr = factor(chr, levels = c(1:22, "X")))

candidate_genes <- coords %>%
    split(.$chr) %>%
    map_dfr(. %>% 
	    inner_join(gencode, by = "chr") %>%
	    filter(tss >= start, tss <= end) %>%
	    distinct(gene_id, gene_name))

write_tsv(candidate_genes, "./sle_candidate_genes.tsv")
