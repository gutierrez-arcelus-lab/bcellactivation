# ==============================================================================
# Description:  Parses the raw QC metadata into a format compatible with the 
#               Salmon Slurm array script. Specifically, it merges technical 
#               replicates/lanes by concatenating FASTQ paths with commas.
# Input:        metadata_longformat.tsv
# Output:       metadata.tsv (Headerless table for Slurm array parsing)
#               Format: Three columns (sample_id, fastq_1, fastq_2) where the 
#               fastq columns are comma-separated lists of all fastqs for a sample.
# ==============================================================================

library(tidyverse)

meta <- 
    "../0_qc/data/metadata_longformat.tsv" |>
    read_tsv() |>
    mutate(sample_id = str_remove(sample_id, "\\.rep\\.\\d$")) |>
    group_by(sample_id, stim, time) |>
    summarise_at(vars(fq1, fq2), 
		 ~paste(., collapse = ",")) |>
    ungroup() |>
    filter(sample_id != "BLANK") |>
    unite("id", c(sample_id, stim, time), sep = "_")

samples_pass <-
    read_tsv("../0_qc/data/samples_pass.tsv")

meta_filtered <- 
    meta |>
    filter(id %in% samples_pass$sample_name)

write_tsv(meta_filtered, "./data/metadata.tsv", col_names = FALSE)
