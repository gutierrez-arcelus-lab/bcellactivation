library(tidyverse)
library(readxl)

patient_data <- 
    "./mgb_data/1. sle_biobank_nodate_n536-deidentified.xlsx" |>
    read_excel() |>
    select(subject_id = Subject_Id, gender = Gender, race = Race1)

control_data <- 
    "./mgb_data/2. controls_nodate_deidentified.xlsx" |>
    read_excel() |>
    select(subject_id = Subject_Id, gender = Gender, race = Race1)

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

sle_data |> count(group, sort = TRUE)

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

# Set up Job array to filter MGB VCFs
dir.create("./mgb_data/ids_per_batch")

expand_grid(batch = batches, chr = paste0("chr", c(1:22, "X"))) |>
    write_tsv("array.spec", col_names = FALSE)

mgb_sle_data |>
    summarise(id = list(sample_id), .by = batch) |>
    mutate(f = paste0("mgb_data/ids_per_batch/", batch, ".txt")) |>
    {function(dat) walk2(dat$id, dat$f, ~write_lines(.x, .y))}()

