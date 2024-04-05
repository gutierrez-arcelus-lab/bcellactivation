library(tidyverse)

assay_dir <- "/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/B_cells_rnaseq"
fq_dir <- file.path(assay_dir, "/34.198.31.178/221025_MG10430_fastq")

meta <- file.path(assay_dir, "metadata.tsv") |>
    read_tsv(col_types = c(.default = "c"))

ids_df <- "../0-genotypes/data/metadata.tsv" |>
    read_tsv(col_types = c(.default = "c"))

# need to know the sample id in the VCF
meta |>
    left_join(ids_df, by = "donor_id") |>
    select(donor_id, sample_id, stim, mgb_id) |>
    write_tsv("./array_spec.tsv", col_names = FALSE)

