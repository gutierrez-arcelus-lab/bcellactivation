library(tidyverse)

summ_stats <- 
    "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Khunsriraksakul/Khunsriraksakul_summstats_mamt.txt" |>
    read_tsv()

pt_vars <- 
    "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Khunsriraksakul/P+T-5e-4_mtag_meta_new.txt" |>
    read_tsv(col_names = c("variant_id", "effect_allele", "beta")) |>
    left_join(select(summ_stats, variant_id, rsid),
	      join_by(variant_id))

dir.create("data")

pt_vars |>
    pull(rsid) |>
    write_lines("./data/pt_5e-4_rsids.txt")


pt_vars |>
    left_join(summ_stats, join_by(variant_id, rsid)) |>
    select(chr, pos = `pos(hg19)`, rsid) |>
    mutate(start = pos - 2.5e5L, stop = pos + 2.5e5L,
	   region = sprintf("%s:%d-%d", chr, start, stop)) |>
    arrange(chr, pos) |>
    select(rsid, region) |>
    write_tsv("./data/pt_5e-4_regions.tsv", col_names = FALSE)
