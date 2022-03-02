library(tidyverse)

bed <- read_tsv("./geuvadis_salmon_quants_ebv.bed")

geuvadis_ids <- read_tsv("./geuvadis_covariates.tsv") %>%
    select(sampleid, kgp_id)

inds <- read_tsv("./qtltools_cov.txt") %>%
    select(-id) %>%
    names()

ebv_bed <- filter(bed, `#chr` == "ebv") %>%
    pivot_longer(starts_with("ERR188"), names_to = "sampleid") %>%
    inner_join(geuvadis_ids) %>%
    select(-sampleid) %>%
    arrange(kgp_id) %>%
    pivot_wider(names_from = kgp_id, values_from = value) %>%
    select(`#chr`, start, end, id, gid, strd, all_of(inds)) %>%
    arrange(start, end)

write_tsv(ebv_bed, "./ebv_phenotypes.bed")
