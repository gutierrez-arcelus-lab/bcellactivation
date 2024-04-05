library(tidyverse)

annotations <- 
    file.path("/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens",
	      "gencode.v39.primary_assembly.annotation.gtf") |>
    read_tsv(comment = "#", col_names = FALSE, col_types = "ccciicccc")

gene_df <- annotations |>
    filter(X3 == "gene", X1 %in% paste0("chr", c(1:22, "X"))) |>
    mutate(gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
	   gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+")) |>
    select(chr = X1, start = X4, end = X5, gene_id, gene_name)

ase_df <- read_tsv("./quasar_results.tsv")

# Some variants are in 2 genes
annot_ase_vars <- 
    ase_df |>
    distinct(snp_id) |>
    separate(snp_id, c("chr", "pos", "ref", "alt"), 
	     sep = ":", convert = TRUE, remove = FALSE) |>
    left_join(gene_df, join_by(chr, between(pos, start, end))) |>
    select(snp_id, gene_id, gene_name)

annot_ase_vars |>
    distinct(snp_id, gene_id, gene_name) |>
    add_count(snp_id) |>
    arrange(desc(n)) |>
    filter(n != 22)

ase_annot_df <- ase_df |>
    left_join(annot_ase_vars, join_by(snp_id), relationship = "many-to-many") |>
    select(donor_id, sample_id, snp_id, gene_id, gene_name, beta:qval)

ase_annot_df |>
    filter(snp_id == "chr11:35208126:T:C")

quasar_genotypes <- read_tsv("./quasar_genotypes.tsv")

quasar_genotypes |>
    filter(snp_id == "chr11:35208126:T:C")
    



