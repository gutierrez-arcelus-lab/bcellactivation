library(tidyverse)

if (!file.exists("data")) dir.create("data")

# Langefeld et al variants
gwas_vars <- 
    "/lab-share/IM-Gutierrez-e2/Public/vitor/lupus/sle_variants/paper_data/langefeld_tier1_hg38.tsv" |>
    read_tsv() |>
    mutate(chr = factor(chr, levels = paste0("chr", 1:22))) |>
    arrange(chr, pos)

write_tsv(gwas_vars, "./data/langefeld_hits.tsv")

gwas_vars |>
    select(chr, pos) |>
    write_tsv("./data/gwas_vars.tsv", col_names = FALSE)

# MGBB array specification
"../../0-genotypes/array_spec.txt" |>
    read_tsv(col_names = FALSE) |>
    filter(X2 != "chrX") |>
    write_tsv("./array_spec.txt", col_names = FALSE)


