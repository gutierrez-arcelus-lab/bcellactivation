library(tidyverse)

batches <- 
    "/lab-share/IM-Gutierrez-e2/Public/External_datasets/mgb_biobank/" |>
    list.files(pattern = "^04")

chrs <- paste0("chr", c(1:22, "X"))

# array for subsetting MGB VCFs
expand_grid(batches, chrs) |>
    write_tsv("./array.spec", col_names = FALSE)


# array for running MBV

read_tsv("../arrayspec_ase.tsv", col_names = FALSE) |>
    select(sample_id = X2, stim = X3) |>
    filter(sample_id %in% c("10044277.1", "10085290.1", "10098629.1")) |>
    unite("bam_prefix", c(sample_id, stim), sep = "_") |>
    expand_grid(batch = batches) |>
    write_tsv("./mbv.spec", col_names = FALSE)

