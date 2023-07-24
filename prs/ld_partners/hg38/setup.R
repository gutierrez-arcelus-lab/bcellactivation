library(tidyverse)
library(readxl)

kgp_pops <- 
    "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/20131219.populations.tsv" |>
    read_tsv() |>
    janitor::clean_names() |>
    filter(population_description != "Total")

kgp_samples <- 
    "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index" |>
    read_tsv(comment = "##") |>
    select(sample_name = SAMPLE_NAME, population_code = POPULATION) |>
    distinct() |>
    left_join(kgp_pops, by = "population_code") |>
    select(sample_name, population_code, super_population)

gwas_ratios <- 
    "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Khunsriraksakul/Khunsriraksakul_SuppData1.xlsx" |>
    read_excel(skip = 2) |>
    select(-Acronym) |>
    fill(Trait, PMID) |>
    mutate(Ancestry = sub("-\\d$", "", Ancestry)) |>
    group_by(Ancestry) |>
    summarise(Neff = sum(Neff)) |>
    ungroup() |>
    filter(Ancestry %in% kgp_pops$super_population) |>
    mutate(p = Neff/sum(Neff)) |>
    arrange(desc(p)) |>
    select(Ancestry, p) |>
    deframe()

kgp_filtered <- kgp_samples |>
    filter((super_population == "EUR") | 
	   (super_population == "EAS" & population_code %in% c("CHB", "JPT")) |
	   (super_population == "AMR" & population_code == "PEL") |
	   (super_population == "SAS")) |>
    mutate(super_population = factor(super_population, levels = names(gwas_ratios)))

max_group_size <- kgp_filtered |> 
    filter(super_population == names(gwas_ratios[which.max(gwas_ratios)])) |>
    nrow()

group_sizes <- round(gwas_ratios / max(gwas_ratios) * max_group_size)

kgp_out <- kgp_filtered |>
    group_split(super_population) |>
    map2_dfr(group_sizes, ~slice_sample(.x, n = .y))

kgp_out |>
    arrange(sample_name) |>
    pull(sample_name) |>
    write_lines("./data/kgp_samples.txt")

# GWAS windows
chr_ids <- 
    "/lab-share/IM-Gutierrez-e2/Public/References/dbSNP/GCF_000001405.39_GRCh38.p13_assembly_report.txt" |>
    read_tsv(comment = "#", col_names = FALSE) |>
    filter(X2 == "assembled-molecule", X8 == "Primary Assembly") |>
    select(nc_name = X7, chr = X10)

rsids <- 
    tibble(rsid = read_lines("/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Khunsriraksakul/rsids.txt"))

gwas_rsids_hg38 <-  
    "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Khunsriraksakul/rsid_dbsnp_hg38.vcf" |>
    read_tsv(comment = "##") |>
    inner_join(chr_ids, join_by(`#CHROM` == nc_name)) |>
    select(chr, POS, ID, REF, ALT) |>
    mutate(chr = factor(chr, levels = sprintf("chr%s", c(1:22, "X", "Y"))))

gwas_rsids_hg38 |>
    select(chr, POS) |>
    write_tsv("./data/gwas_positions.tsv", col_names = FALSE)

#unmapped <- anti_join(rsids, gwas_rsids_hg38, join_by(rsid == ID))

sentinel <- 
    "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Khunsriraksakul/Khunsriraksakul_SuppData2.xlsx" |>
    read_excel(skip = 2) |>
    janitor::clean_names() |>
    left_join(gwas_rsids_hg38, join_by(rs_id == ID)) |>
    select(chr, pos = POS, rs_id, effect_allele, other_allele, 
	   REF, ALT, beta, se, p_value, mapped_gene) |>
    arrange(chr, pos)

extra_sentinel <- 
    tibble(rs_id = read_lines("../data/ExtraSentinels_500k.tsv")) |>
    inner_join(gwas_rsids_hg38, join_by(rs_id == ID)) |>
    select(chr, pos = POS, rs_id)

all_sentinels <- 
    select(sentinel, chr, pos, rs_id) |>
    bind_rows(extra_sentinel) |>
    arrange(chr, pos) |>
    mutate(start = pos - 5e5L, stop = pos + 5e5L,
	   region = sprintf("%s:%d-%d", chr, start, stop))

all_sentinels |>
    select(region, rs_id) |>
    write_tsv("./data/regions.tsv", col_names = FALSE)

# summary stats
summ_stats_hg19 <- 
    "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Khunsriraksakul/Khunsriraksakul_summstats_mamt.txt" |>
    read_tsv() |>
    janitor::clean_names()

summ_stats_hg38 <- 
    summ_stats_hg19 |>
    select(rsid, effect_allele, other_allele, freq, beta, se, z, pval) |>
    inner_join(gwas_rsids_hg38, join_by(rsid == ID)) |>
    select(chr, pos = POS, rsid, ref = REF, alt = ALT, effect_allele, other_allele, freq:pval) |>
    arrange(chr, pos)

write_tsv(summ_stats_hg38, "./data/Khunsriraksakul_summstats_hg38.tsv")
