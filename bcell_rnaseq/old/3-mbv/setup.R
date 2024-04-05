library(tidyverse)

batches <- 
    "/lab-share/IM-Gutierrez-e2/Public/External_datasets/mgb_biobank/" |>
    list.files(pattern = "^04")

chrs <- paste0("chr", c(1:22, "X"))

# array for subsetting MGB VCFs
array_file <- "array_spec_vcf.txt"
array_mbv <- "array_spec_mbv.txt"

expand_grid(batches, chrs) |>
    write_tsv(array_file, col_names = FALSE)

# array for running MBV
read_tsv("../2-ase/array_spec.tsv", col_names = FALSE) |>
    select(sample_id = X2, stim = X3) |>
    unite("bam_prefix", c(sample_id, stim), sep = "_") |>
    expand_grid(batch = batches) |>
    write_tsv(array_mbv, col_names = FALSE)

