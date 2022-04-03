library(tidyverse)
library(rhdf5)

# SNP LOC
gwas <- read_tsv("./data/bentham_2015_26502338_sle_efo0002690_1_gwas.sumstats.tsv.gz")

snp_loc <- select(gwas, rsid, chrom, pos)
pval <- select(gwas, SNP = rsid, P = p)


write_tsv(snp_loc, "./data/bentham_snp_loc.tsv", col_names = FALSE)
write_tsv(pval, "./data/bentham_pval.tsv", col_names = TRUE)

# GENE LOC
sc_genes <- h5read("./data/bcells_singlet.h5ad", "var/gene")

gene_strand <- read_tsv("../../data/gencodev38_genestrand.tsv")

gencode <- read_tsv("../../data/gencodev38_genes_hg19.bed", col_names = FALSE) %>%
    extract(X4, c("gene_id", "gene_name"), "(ENSG[0-9.]+)_(.+)") %>%
    mutate(chr = sub("^chr", "", X1)) %>%
    left_join(gene_strand, by = "gene_id") %>%
    select(gene_name, chr, start = X2, end = X3, strand, gene_id) %>%
    mutate(gene_id = sub("\\.\\d+$", "", gene_id)) %>%
    filter(gene_id %in% sc_genes) %>%
    rownames_to_column("i") %>%
    unite(gene_name, c("i", "gene_name"), sep = "_") %>%
    select(gene_id, chr, start, end, strand, gene_name)


write_tsv(gencode, "./data/gencode_genes/genes.v38.hg19.loc", col_names = FALSE)
