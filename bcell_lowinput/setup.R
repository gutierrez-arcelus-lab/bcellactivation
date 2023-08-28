library(tidyverse)

meta <- 
    "/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/B_cells_lowinput/metadata/description_hsapiens_qc.tsv" |>
    read_tsv() |>
    mutate(sample_id = str_remove(sample_id, "\\.rep\\.\\d$")) |>
    group_by(sample_id, stim, time) |>
    summarise_at(vars(fq1, fq2), ~paste(., collapse = ",")) |>
    ungroup() |>
    filter(sample_id != "BLANK") |>
    unite("id", c(sample_id, stim, time), sep = "_")

write_tsv(meta, "./data/metadata_pooledreps.tsv", col_names = FALSE)
