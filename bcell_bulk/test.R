library(tidyverse)


conditions <- c("16hr_resting", "24hr_IgG", "72hr_IgG", "24hr_RSQ", "72hr_RSQ")

ase <- read_tsv("./results/phaser/MG8989_gene_ae.txt") %>%
    mutate(bam = sub("20210615_(\\d+hr_[^_]+)_.+$", "\\1", bam),
	   bam = factor(bam, levels = conditions)) %>%
    extract(name, c("gene_id", "gene_name"), "(ENSG[0-9.]+)_(.+)") %>%
    select(condition = bam, everything()) %>%
    arrange(contig, start, gene_name, condition)

ase %>%
    mutate(ai = abs(0.5 - aCount/totalCount)) %>%
    group_by(gene_name) %>%
    filter(all(!is.na(variants))) %>%
    filter(any(totalCount >= 50)) %>%
    group_by(gene_name, variants) %>%
    filter(first(ai) < .1 & any(ai > .2 & totalCount >= 16)) %>%
    mutate(score = ai - first(ai)) %>%
    summarise(score = max(score)) %>%
    separate_rows(variants, sep = ",") %>%
    filter(n() >= 5) %>%
    summarise(variants = paste(variants, collapse = ","),
	      score = unique(score)) %>%
    ungroup() %>%
    arrange(desc(score)) %>%
    print(n = 50)

ase %>%
    filter(gene_name == "C22orf34") %>%
    mutate(ai = abs(0.5 - aCount/totalCount)) %>%
    mutate(score = ai - first(ai))

