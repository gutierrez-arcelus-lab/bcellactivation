library(tidyverse)

# GRCh38
annotations <- "../data/gencode.v38.primary_assembly.annotation.gtf" %>%
    read_tsv(comment = "#", col_names = FALSE, col_types = "ccciicccc") %>%
    filter(X3 == "gene", X1 %in% paste0("chr", c(1:22, "X"))) %>%
    mutate(gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
	   gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+"),
	   gene_type = str_extract(X9, "(?<=gene_type\\s\")[^\"]+")) %>%
    select(chr = 1, start = 4, end = 5, gene_id, gene_name, gene_type, strand = X7)

annotations %>%
    select(-strand, -gene_type) %>%
    unite("gene", c("gene_id", "gene_name"), sep = "_") %>%
    write_tsv("./gencodev38_genes.bed", col_names = FALSE)

annotations %>%
    select(gene_id, strand) %>%
    write_tsv("./gencodev38_genestrand.tsv")

annotations %>%
    select(gene_id, gene_name, gene_type) %>%
    write_tsv("./gencode_v38_genetypes.tsv")

# hg19

annot_hg19 <- "../data/gencode.v19.chr_patch_hapl_scaff.annotation.gtf" %>%
    read_tsv(comment = "#", col_names = FALSE, col_types = "ccciicccc") %>%
    filter(X3 == "gene", X1 %in% paste0("chr", c(1:22, "X"))) %>%
    mutate(gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
	   gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+")) %>%
    select(chr = 1, start = 4, end = 5, gene_id, gene_name, strand = X7)

write_tsv(annot_hg19, "./gencode.v19.genes.tsv")
