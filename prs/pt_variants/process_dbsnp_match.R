library(tidyverse)

chr_ids <- 
    "/lab-share/IM-Gutierrez-e2/Public/References/dbSNP/GCF_000001405.39_GRCh38.p13_assembly_report.txt" |>
    read_tsv(comment = "#", col_names = FALSE) |>
    filter(X2 == "assembled-molecule", X8 == "Primary Assembly") |>
    select(nc_name = X7, chr = X10)

vcf <- 
    "./data/pt_5e-4.vcf" |>
    read_tsv(comment = "##") |>
    inner_join(chr_ids, join_by(`#CHROM` == nc_name)) |>
    select(chr, pos = POS, rsid = ID, REF, ALT)
    
vcf |>
    select(chr, pos) |>
    write_tsv("./data/pt_5e-4_grch38_positions.tsv", col_names = FALSE)

vcf |>
    write_tsv("./data/pt_5e-4_grch38.tsv")
