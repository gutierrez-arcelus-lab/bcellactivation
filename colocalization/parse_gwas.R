library(tidyverse)

bentham_stats <- "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Bentham/ebi-a-GCST003156.vcf.gz" |>
    read_tsv(comment = "##")

bentham_region <- bentham_stats |>
    filter(`#CHROM` == 1, between(POS, 205e6, 208e6)) |>
    separate(`EBI-a-GCST003156`, c("eff", "se", "logp", "id"), sep = ":") |>
    select(id = ID, pos = POS, logp)

langefeld_stats_chr1 <- "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Langefeld/ea.imputed.chr1.out" |>
    read_delim(comment = "#", delim = " ")

langefeld_stats_chr7 <- "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Langefeld/ea.imputed.chr7.out" |>
    read_delim(comment = "#", delim = " ")

langefeld_region <- langefeld_stats_chr1 |>
    filter(between(position, 205e6, 208e6)) |>
    arrange(frequentist_add_wald_pvalue_1) |> 
    select(1:2, 4, frequentist_add_wald_pvalue_1, cases_maf, controls_maf)

write_tsv(bentham_region, "./data/bentham_chr1.tsv")
write_tsv(langefeld_region, "./data/langefeld_chr1.tsv")

lange_vars_chr1 <- langefeld_stats_chr1 |>
    filter(grepl("^rs\\d+", rsid)) |>
    extract(rsid, c("rsid"), "(^rs\\d+)") |>
    distinct(rsid) |>
    pull(rsid)
    
lange_vars_chr7 <- langefeld_stats_chr7 |>
    filter(grepl("^rs\\d+", rsid)) |>
    extract(rsid, c("rsid"), "(^rs\\d+)") |>
    distinct(rsid) |>
    pull(rsid)

bentham_vars <- bentham_stats |>
    filter(`#CHROM` %in% c(1, 7)) |>
    distinct(rsid = ID) |>
    pull(rsid)

all_variants <- unique(c(lange_vars_chr1, lange_vars_chr7, bentham_vars))

write_lines(sort(all_variants), "./data/chr1_chr7_EURGWAS_variants.txt")

