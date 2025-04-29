library(tidyverse)

meta <- 
    "../0-qc_rnaseq/metadata.tsv" |>
    read_tsv(col_types = "ccccc") |>
    separate_rows(R1:R2, sep = ",")

write_tsv(meta, "./metadata.tsv")

# MBV
mbv <- 
    "../2-mbv/matching_results.tsv" |>
    read_tsv(col_types = "ccccc") |>
    separate(rnaseq_sample_id, c("donor_id", "replic_id"), sep = "_") |>
    select(donor_id, replic_id, stim = rnaseq_stim, vcf_donor_id)

meta_qced <- 
    meta |>
    mutate_at(vars(R1:R2), basename) |>
    mutate(R1 = str_replace(R1, "\\.fastq\\.gz", "_val_1.fq.gz"),
	   R2 = str_replace(R2, "\\.fastq\\.gz", "_val_2.fq.gz")) |>
    mutate_at(vars(R1:R2), ~file.path("/temp_work/ch229163/fastq/highinput", .)) |>
    group_by(donor_id, replic_id, stim) |>
    summarise_at(vars(R1:R2), 
		 ~paste(., collapse = ",")) |>
    ungroup() |>
    left_join(mbv, join_by(donor_id, replic_id, stim)) |>
    unite("sample_id", c(stim, vcf_donor_id, replic_id), sep = "_") |>
    select(sample_id, R1, R2)

write_tsv(meta_qced, "./metadata_qced.tsv")
