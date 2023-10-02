library(tidyverse)
library(readxl)

dir.create("data")

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
summ_stats <- 
    "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Khunsriraksakul/Khunsriraksakul_summstats_mamt.txt" |>
    read_tsv() |>
    janitor::clean_names() |>
    select(chr, pos = pos_hg19, rs_id = rsid)

sentinel <- 
    "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Khunsriraksakul/Khunsriraksakul_SuppData2.xlsx" |>
    read_excel(skip = 2) |>
    janitor::clean_names() |>
    separate(chr_pos_hg19, c("chr", "pos"), sep = ":", convert = TRUE) |>
    select(chr, pos, rs_id)

extra_sentinel <- 
    tibble(rs_id = read_lines("../data/ExtraSentinels_500k.tsv")) |>
    left_join(summ_stats, join_by(rs_id)) |>
    select(chr, pos, rs_id) 

all_sentinels <- 
    bind_rows(sentinel, extra_sentinel) |>
    mutate(chr = factor(chr, levels = c(1:22, "X", "Y"))) |>
    arrange(chr, pos) |>
    mutate(start = pos - 5e5L, stop = pos + 5e5L,
	   region = sprintf("%s:%d-%d", chr, start, stop))

all_sentinels |>
    select(region, rs_id) |>
    write_tsv("./data/regions.tsv", col_names = FALSE)

summ_stats |>
    select(chr, pos) |>
    arrange(chr, pos) |>
    write_tsv("./data/gwas_positions.tsv", col_names = FALSE)


# Only original sentinels






