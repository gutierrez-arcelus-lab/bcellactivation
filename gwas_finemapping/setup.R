library(tidyverse)
library(readxl)
library(glue)

if (!file.exists("data")) dir.create("data")

# TGP
"http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel" |>
    read_table() |>
    filter(super_pop == "EUR") |>
    pull(sample) |>
    write_lines("./data/tgp_eur.txt")

# Langefeld summary stats
parse_pvalues <- function(p) {

    p |>
        {function(x) gsub("\\s", "", x)}() |>
        {function(x) sub("^([0-9.]+)[^0-9.]10[^0-9.](\\d+)$", "\\1e-\\2", x)}() |>
        parse_number()
}

langefeld_ea <-
    "../2-liftover/data/langefeld_tableS2.xlxs" |>
    read_excel(1, skip = 2) |>
    select(rsid = 1,
           chr = 2,
           pos = 3,
           gene_region = 4,
           region_rank = 5,
           ref_allele = 6,
           pval = `P-value`) |>
    filter(!is.na(pos)) |>
    mutate(rsid = gsub("\\s|[a-zA-Z,]+$", "", rsid),
           chr = sub("^(\\d+)[qp]\\d+$", "\\1", chr),
           chr = factor(chr, levels = c(1:22, "X", "Y")),
           pval = parse_pvalues(pval)) |>
    arrange(chr, pos)

target_region <- 
    langefeld_ea |>
    filter(grepl("TNFSF4", gene_region)) |>
    slice_min(pval, n = 1) |>
    mutate(START = pos - 1e6L, END = pos + 1e6L) |>
    select(chr, START, END)

glue("{target_region$chr}:{target_region$START}-{target_region$END}") |>
    write_lines("./data/target_region.txt")

langefeld_summ_stats <-
    "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Langefeld/ea.imputed.chr{target_region$chr}.out" |>
    glue() |>
    data.table::fread(skip = "alternate_ids") |>
    as_tibble() |>
    add_column(chr = target_region$chr, .before = 1) |>
    select(chr,
	   rsid,
	   pos = position,
	   alleleA,
	   alleleB,
	   beta = `frequentist_add_beta_1:add/sle=1`,
	   se = frequentist_add_se_1,
	   p = frequentist_add_lrt_pvalue,
	   ) |>
    mutate(rsid = str_extract(rsid, "rs\\d+")) |>
    inner_join(target_region, join_by(chr, between(pos, START, END))) |>
    select(-START, -END)

langefeld_munge <- 
    langefeld_summ_stats |>
    drop_na(p, beta, se) |>
    mutate(chr = paste0("chr", chr),
           chr = fct_inorder(chr),
	   rsid = replace_na(rsid, ".")) |>
    select(CHR = chr, BP = pos, SNP = rsid, A1 = alleleA, A2 = alleleB, BETA = beta, SE = se, P = p)

write_tsv(langefeld_munge, "./data/langefeld_munge.tsv")

tibble(c1 = names(langefeld_munge), c2 = names(langefeld_munge)) |>
    write_tsv("./data/munge_col_header.tsv", col_names = FALSE)
