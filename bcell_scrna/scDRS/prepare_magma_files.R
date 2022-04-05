library(tidyverse)
library(rhdf5)

# SNP LOC
gwas <- 
    "./data/summ_stats/bentham_2015_26502338_sle_efo0002690_1_gwas.sumstats.tsv.gz" %>%
    read_tsv() %>%
    filter(!(chrom == 6 & between(pos, 29e6, 34e6))) 

snp_loc <- select(gwas, rsid, chrom, pos)
pval <- select(gwas, SNP = rsid, P = p)

write_tsv(snp_loc, "./data/summ_stats/bentham_snp_loc.tsv", col_names = FALSE)
write_tsv(pval, "./data/summ_stats/bentham_pval.tsv", col_names = TRUE)

# GENE LOC
sc_genes <- h5read("./data/expression/bcells_singlet_seurat.h5ad", "/var/_index")

gencode <- read_tsv("../../data/gencode.v19.genes.tsv") %>%
    mutate(chr = sub("^chr", "", chr),
	   chr = recode(chr, X = "23")) %>%
    select(gene_name, chr, start, end, strand, gene_id) %>%
    mutate(gene_id = sub("\\.\\d+$", "", gene_id)) %>%
    filter(gene_id %in% sc_genes) %>%
    rownames_to_column("i") %>%
    unite(gene_name, c("i", "gene_name"), sep = "_") %>%
    select(gene_id, chr, start, end, strand, gene_name)

write_tsv(gencode, "./data/gencode_genes/genes.v19.hg19.loc", col_names = FALSE)
