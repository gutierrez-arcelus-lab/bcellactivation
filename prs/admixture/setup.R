library(tidyverse)

#dir.create("data")

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

kgp_samples |>
    arrange(sample_name) |>
    mutate(id = paste(sample_name, sample_name, sep = "_")) |>
    pull(id) |>
    write_lines("./data/kgp_ids.txt")

kgp_samples |>
    arrange(super_population, population_code, sample_name) |>
    write_tsv("./data/kgp_metadata.tsv")


mgb_ids <-
    list.files("../mgb_data/ids_per_batch", full.names = TRUE) |>
    map(read_lines) |>
    unlist() |>
    sort()

write_lines(mgb_ids, "./data/mgb_ids.txt")
