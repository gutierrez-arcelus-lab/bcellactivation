library(tidyverse)

chr_ids <- 
    "/lab-share/IM-Gutierrez-e2/Public/References/dbSNP/GCF_000001405.39_GRCh38.p13_assembly_report.txt" |>
    read_tsv(comment = "#", col_names = FALSE) |>
    filter(X2 == "assembled-molecule", X8 == "Primary Assembly") |>
    select(nc_name = X7, chr = X10)

sentinel_info <- 
    "./data/Khunsriraksakul_sentinel_r2.tsv" |>
    read_tsv() |>
    select(-r2) |>
    mutate(sentinel = fct_inorder(sentinel)) |>
    group_by(sentinel) |>
    mutate(ld_group = cur_group_id()) |>
    ungroup() |>
    mutate(sentinel = as.character(sentinel)) |>
    select(sentinel, ld_partner = snp2, ld_group) |>
    pivot_longer(sentinel:ld_partner, values_to = "rsid", names_to = "type") |>
    drop_na() |>
    distinct()

sentinels_vcf <- 
    "./data/sentinels_dbsnp.vcf" |>
    read_tsv(comment = "##") |>
    inner_join(chr_ids, join_by(`#CHROM` == nc_name)) |>
    select(chr, pos = POS, rsid = ID, REF, ALT) |>
    left_join(sentinel_info, join_by(rsid))

write_tsv(sentinels_vcf, "./data/sentinels_hg38_positions.tsv")

sentinels_vcf |>
    select(chr, pos) |>
    write_tsv("./data/sentinels_grch38.tsv", col_names = FALSE)
