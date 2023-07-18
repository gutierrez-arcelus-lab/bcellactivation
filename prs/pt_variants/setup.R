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
