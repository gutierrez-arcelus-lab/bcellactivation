library(tidyverse)

annotations <- "../data/gencode.v38.primary_assembly.annotation.gtf" %>%
    read_tsv(comment = "#", col_names = FALSE, col_types = "ccciicccc") %>%
    filter(X3 == "gene", X1 %in% paste0("chr", c(1:22, "X"))) %>%
    mutate(gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
	   gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+")) %>%
    select(chr = 1, start = 4, end = 5, gene_id, gene_name)

annotations %>%
    unite("gene", c("gene_id", "gene_name"), sep = "_") %>%
    write_tsv("./gencodev38_genes.bed", col_names = FALSE)
