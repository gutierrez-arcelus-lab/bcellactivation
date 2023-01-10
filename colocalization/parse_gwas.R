library(tidyverse)

bentham_stats <- "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Bentham/ebi-a-GCST003156.vcf.gz" |>
    read_tsv(comment = "##") |>

bentham_region <- bentham_stats |>
    filter(`#CHROM` == 1, between(POS, 205e6, 208e6)) |>
    separate(`EBI-a-GCST003156`, c("eff", "se", "logp", "id"), sep = ":") |>
    select(id = ID, pos = POS, logp)

langefeld_stats <- "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Langefeld/ea.imputed.chr1.out" |>
    read_delim(comment = "#", delim = " ")

langefeld_region <- langefeld_stats |>
    filter(between(position, 205e6, 208e6)) |>
    arrange(frequentist_add_wald_pvalue_1) |> 
    select(1:2, 4, frequentist_add_wald_pvalue_1, cases_maf, controls_maf)

write_tsv(bentham_region, "./data/bentham_chr1.tsv")
write_tsv(langefeld_region, "./data/langefeld_chr1.tsv")
