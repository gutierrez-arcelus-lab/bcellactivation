library(tidyverse)

# Create output directory
if (!file.exists("data")) dir.create("data")

# RNA-seq metadata
assay_dir <- "/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/B_cells_rnaseq"

meta <- file.path(assay_dir, "metadata.tsv") |>
    read_tsv(col_types = c(.default = "c"))

# MGB metadata
batches <- sprintf("04%02d", 1:8)

ids_df <- 
    "/lab-share/IM-Gutierrez-e2/Public/External_datasets/mgb_biobank/%s/etc/ids.tsv" |>
    sprintf(batches) |>
    setNames(batches) |>
    map_df(~read_tsv(., col_names = c("donor_id", "mgb_id"), col_types = "c"), .id = "batch") |>
    filter(donor_id %in% unique(meta$donor_id))

write_tsv(ids_df, "./data/metadata.tsv")

ids_df |>
    select(batch, mgb_id) |>
    group_by(batch) |>
    summarise(id = list(mgb_id)) |>
    ungroup() |>
    mutate(f = paste0("data/", batch, ".txt")) |>
    {function(dat) walk2(dat$id, dat$f, ~write_lines(.x, .y))}()

# Slurm array spec for VCF processing
ids_df |>
    distinct(batch) |>
    cross_join(tibble(chr = sprintf("chr%s", c(1:22, "X")))) |>
    write_tsv("./array_spec.txt", col_names = FALSE)

