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
    mutate_at(vars(R1:R2), ~file.path("/temp_work/ch229163/fastq/highinput", .)) |>
    group_by(donor_id, replic_id, stim) |>
    summarise_at(vars(R1:R2), 
		 ~paste(., collapse = ",")) |>
    ungroup() |>
    left_join(mbv, join_by(donor_id, replic_id, stim)) |>
    unite("sample_id", c(stim, vcf_donor_id, replic_id), sep = "_") |>
    select(sample_id, R1, R2)

write_tsv(meta_qced, "./data/metadata_qced.tsv")

# Leafcutter
# remove sample 10028815_unstday0 because it has a low proportion of uniquely mapped reads
# if not removed, it is an outlier in the PCA by leafviz

meta_qced_filt <-
    meta_qced |>
    select(sample_id) |>
    filter(sample_id != "unstday0_10028815_1")

meta_qced_filt |>
    mutate(juncfile = paste0(getwd(), "/bam/", sample_id, ".junc")) |>
    pull(juncfile) |>
    write_lines("./data/junction_files.txt")

# make groups files
samples_for_ds <- 
    meta_qced_filt |>
    separate(sample_id, c("stim", "donor", "replic"), sep = "_", remove = FALSE) |>
    filter(replic == 1) 

contrasts_df <- 
    combn(c("unstday0", "TLR7", "BCR", "DN2"), 2) |>
    t() |>
    as.data.frame() |>
    mutate(contrast = paste0(V1, "vs.", V2)) |>
    select(contrast, stim1 = V1, stim2 = V2) |>
    pivot_longer(-contrast) |>
    select(-name)

left_join(contrasts_df, samples_for_ds, 
	  join_by(value == stim), 
	  relationship = "many-to-many") |>
    select(contrast, sample_id, group = value) |>
    group_by(contrast) |>
    nest() |>
    ungroup() |>
    mutate(file_name = glue::glue("./data/groups_{contrast}.tsv")) |>
    select(data, file_name) |>
    pwalk(~write_tsv(.x, .y, col_names = FALSE))

# save contrasts for slurm
contrasts_df |>
    distinct(contrast) |>
    pull(contrast) |>
    write_lines("./data/contrasts.txt")

