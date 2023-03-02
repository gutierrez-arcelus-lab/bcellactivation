library(tidyverse)

samples <- read_lines("./data/mgb_sample_ids.txt")

batches <- sprintf("04%02d", 1:8)

ids_df <- 
    "/lab-share/IM-Gutierrez-e2/Public/External_datasets/mgb_biobank/%s/etc/ids.tsv" |>
    sprintf(batches) |>
    setNames(batches) |>
    map_df(~read_tsv(., col_names = c("subject_id", "mgb_id")), .id = "batch") |>
    filter(subject_id %in% samples)

ids_df |>
    select(batch, mgb_id) |>
    group_by(batch) |>
    summarise(id = list(mgb_id)) |>
    ungroup() |>
    mutate(f = paste0("./data/", batch, ".txt")) |>
    {function(dat) walk2(dat$id, dat$f, ~write_lines(.x, .y))}()

write_tsv(ids_df, "./data/mgb_sample_info.tsv", col_names = FALSE)

