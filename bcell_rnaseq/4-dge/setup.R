library(tidyverse)

if (!file.exists("data")) dir.create("data")
if (!file.exists("results")) dir.create("results")

# Exclude library 231227 because some files might be corrupt
meta <- 
    "../0-qc_rnaseq/metadata.tsv" |>
    read_tsv(col_types = "ccccc") |>
    separate_rows(R1:R2, sep = ",") |>
    filter(!grepl("/231227", R1))

write_tsv(meta, "./data/metadata.tsv")

# MBV
mbv <- 
    "../2-mbv/matching_results.tsv" |>
    read_tsv(col_types = "ccccc") |>
    separate(rnaseq_sample_id, c("donor_id", "replic_id"), sep = "_") |>
    select(donor_id, replic_id, stim = rnaseq_stim, vcf_donor_id)

# QCed fastqs for mapping
meta_qced <- 
    meta |>
    mutate_at(vars(R1:R2), basename) |>
    mutate(R1 = str_replace(R1, "\\.fastq\\.gz", "_val_1.fq.gz"),
	   R2 = str_replace(R2, "\\.fastq\\.gz", "_val_2.fq.gz")) |>
    mutate_at(vars(R1:R2), ~file.path("/temp_work/ch229163/fastq", .)) |>
    group_by(donor_id, replic_id, stim) |>
    summarise_at(vars(R1:R2), 
		 ~paste(., collapse = ",")) |>
    ungroup() |>
    left_join(mbv, join_by(donor_id, replic_id, stim)) |>
    unite("sample_id", c(stim, vcf_donor_id, replic_id), sep = "_") |>
    select(sample_id, R1, R2)

write_tsv(meta_qced, "./data/metadata_qced.tsv")

# Collapse technical replicates for salmon_poolreps analysis
meta_qced_pooled <- 
    meta |>
    mutate_at(vars(R1:R2), basename) |>
    mutate(R1 = str_replace(R1, "\\.fastq\\.gz", "_val_1.fq.gz"),
	   R2 = str_replace(R2, "\\.fastq\\.gz", "_val_2.fq.gz")) |>
    mutate_at(vars(R1:R2), ~file.path("/temp_work/ch229163/fastq", .)) |>
    group_by(donor_id, replic_id, stim) |>
    summarise_at(vars(R1:R2), ~paste(., collapse = ",")) |>
    ungroup() |>
    mutate(stim = recode(stim,
			 "unstday0" = "unst.0", 
			 "BCR" = "BCR.24", 
			 "TLR7" = "TLR7.24",
			 "DN2" = "DN2.24"),
	   stim = factor(stim, levels = c("unst.0", "TLR7.24", "BCR.24", "DN2.24"))) |>
    arrange(donor_id, replic_id, stim) |>
    group_by(donor_id, stim) |>
    summarise_at(vars(R1:R2), ~paste(., collapse = ",")) |>
    ungroup() |>
    unite("sample_id", c(stim, donor_id), sep = ".") |>
    select(sample_id, R1, R2)
    
write_tsv(meta_qced_pooled, "./data/metadata_qced_pooled.tsv")

