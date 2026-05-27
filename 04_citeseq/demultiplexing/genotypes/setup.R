library(tidyverse)

samples <- read_lines("./data/mgb_sample_ids.txt")

batches <- sprintf("04%02d", 1:8)

ids_df <- 
    "/lab-share/IM-Gutierrez-e2/Public/External_datasets/mgb_biobank/{batches}/etc/ids.tsv" |>
    glue::glue() |>
    setNames(batches) |>
    map_df(~read_tsv(., col_names = c("subject_id", "mgb_id")), .id = "batch") |>
    filter(subject_id %in% samples)

ids_df |>
    group_by(batch) |>
    summarise(sample_ids = paste(mgb_id, collapse = ",")) |>
    ungroup() |>
    expand_grid(tibble(chrom = c(1:22, "X"))) |>
    select(chrom, batch, sample_ids) |>
    write_tsv("array_spec.txt", col_names = FALSE)
