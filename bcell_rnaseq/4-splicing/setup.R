library(tidyverse)

# Exclude library 231227 because some files might be corrupt
meta <- 
    "../0-qc_rnaseq/metadata.tsv" |>
    read_tsv(col_types = "ccccc") |>
    separate_rows(R1:R2, sep = ",") |>
    filter(!grepl("/231227", R1))

write_tsv(meta, "./metadata.tsv")

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
    mutate_at(vars(R1:R2), ~file.path("/temp_work/ch229163/fastq/highinput", .)) |>
    group_by(donor_id, replic_id, stim) |>
    summarise_at(vars(R1:R2), 
		 ~paste(., collapse = ",")) |>
    ungroup() |>
    left_join(mbv, join_by(donor_id, replic_id, stim)) |>
    unite("sample_id", c(stim, vcf_donor_id, replic_id), sep = "_") |>
    select(sample_id, R1, R2)

write_tsv(meta_qced, "./metadata_qced.tsv")

# Leafcutter
if (!file.exists("data")) dir.create("data")
if (!file.exists("results")) dir.create("results")

meta_qced |>
    select(sample_id) |>
    mutate(juncfile = paste0(getwd(), "/bam/", sample_id, ".junc")) |>
    pull(juncfile) |>
    write_lines("./data/junction_files.txt")

groups_file_tlr7 <-
    meta_qced |>
    select(sample_id) |>
    separate(sample_id, c("stim", "donor", "replic"), sep = "_", remove = FALSE) |>
    filter(replic == 1) |>
    filter(stim %in% c("unstday0", "TLR7")) |>
    mutate(stim = factor(stim, levels = c("unstday0", "TLR7"))) |>
    arrange(stim, donor, replic) |>
    select(sample_id, stim)

groups_file_bcr <-
    meta_qced |>
    select(sample_id) |>
    separate(sample_id, c("stim", "donor", "replic"), sep = "_", remove = FALSE) |>
    filter(replic == 1) |>
    filter(stim %in% c("unstday0", "BCR")) |>
    mutate(stim = factor(stim, levels = c("unstday0", "BCR"))) |>
    arrange(stim, donor, replic) |>
    select(sample_id, stim)

groups_file_dn2 <-
    meta_qced |>
    select(sample_id) |>
    separate(sample_id, c("stim", "donor", "replic"), sep = "_", remove = FALSE) |>
    filter(replic == 1) |>
    filter(stim %in% c("unstday0", "DN2")) |>
    mutate(stim = factor(stim, levels = c("unstday0", "DN2"))) |>
    arrange(stim, donor, replic) |>
    select(sample_id, stim)

write_tsv(groups_file_tlr7, "./data/groups_TLR7.txt", col_names = FALSE)
write_tsv(groups_file_bcr, "./data/groups_BCR.txt", col_names = FALSE)
write_tsv(groups_file_dn2, "./data/groups_DN2.txt", col_names = FALSE)
