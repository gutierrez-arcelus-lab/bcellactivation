library(tidyverse)

datasets <- 
    "./data/eqtl_datasets.tsv" |>
    read_tsv() |>
    mutate(tissue_label = case_when(condition_label == "naive" ~ tissue_label,
				    TRUE ~ paste(tissue_label, condition_label))) |>
    select(dataset_id, study = study_label, tissue = tissue_label, mol_phenotype = quant_method)

temp_dir <- system("echo $TEMP_WORK", intern = TRUE) |>
    file.path("eqtl_catalogue")

query_df <- list.files(temp_dir, full.names = TRUE, pattern = "chr12") |>
    {function(x) setNames(x, sub("^([^_]+).+$", "\\1", basename(x)))}() |>
    map_dfr(read_tsv, .id = "dataset_id")

left_join(query_df, datasets, join_by(dataset_id)) |>
    select(names(datasets), everything()) |>
    write_tsv("./data/eqtl_catalogue_associations.tsv.gz")
