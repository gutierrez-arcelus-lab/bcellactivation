library(tidyverse)
library(readxl)

patient_data <- 
    "./mgb_data/1. sle_biobank_nodate_n536-deidentified.xlsx" |>
    read_excel() |>
    janitor::clean_names() |>
    select(subject_id, gender, race = race1, ethnic_group)

control_data <- 
    "./mgb_data/2. controls_nodate_deidentified.xlsx" |>
    read_excel() |>
    janitor::clean_names() |>
    select(subject_id, gender, race = race1, ethnic_group)

sle_data <- 
    bind_rows("SLE" = patient_data, 
	      "Control" = control_data, 
	      .id = "group") |>
    mutate(race = case_match(race, 
			     c("Declined", "Not Given", "Not Reported", "Unavailable", "Unknown") ~ NA,
			     "American Indian" ~ "American Indian or Alaska Native",
			     "Hispanic" ~ "Hispanic or Latino",
			     "Hawaiian" ~ "Native Hawaiian or Other Pacific Islander",
			     .default = race))

write_tsv(sle_data, "./mgb_data/sle_data.tsv")


batches <- sprintf("04%02d", c(1:8, 10))

mgb_ids <- 
    "/lab-share/IM-Gutierrez-e2/Public/External_datasets/mgb_biobank/%s/etc/ids.tsv" |>
    sprintf(batches) |>
    setNames(batches) |>
    map_df(~read_tsv(., col_names = c("subject_id", "sample_id")), .id = "batch")

# Some individuals in early batches were genotyped again in 0410
# Keep just first record for each individual
mgb_sle_data <- 
    inner_join(sle_data, mgb_ids, join_by(subject_id)) |>
    add_count(subject_id) |>
    filter(n == 1 | (n == 2 & batch != "0410")) |>
    select(-n)

count(mgb_sle_data, group)
count(mgb_sle_data, batch)
count(mgb_sle_data, race)

# Set up Job array to filter MGB VCFs
#dir.create("./mgb_data/ids_per_batch")

expand_grid(batch = batches, chr = paste0("chr", c(1:22, "X"))) |>
    write_tsv("array.spec", col_names = FALSE)

mgb_sle_data |>
    summarise(id = list(sample_id), .by = batch) |>
    mutate(f = paste0("mgb_data/ids_per_batch/", batch, ".txt")) |>
    {function(dat) walk2(dat$id, dat$f, ~write_lines(.x, .y))}()

# 1000 Genomes unrelated samples
kgp <- 
    "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index" |>
    read_tsv(comment = "##") |>
    pull(SAMPLE_NAME)

write_lines(kgp, "./mgb_data/kgp_samples.txt")

# ADMIXTURE ref panel
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
