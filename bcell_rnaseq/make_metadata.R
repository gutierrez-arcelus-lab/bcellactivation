library(tidyverse)

fqdir <- "/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/B_cells_rnaseq/34.198.31.178/221025_MG10430_fastq/"

meta <- tibble(fq = list.files(fqdir, pattern = "fastq\\.gz$")) |>
    mutate(info = sub("^\\d+_(\\d+[_123]*)_([^_]+)_.+$", "\\1-\\2", fq)) |>
    separate(info, c("subject_id", "stim"), sep = "-") |>
    extract(fq, "r", ".+_(R[12])_.+", remove = FALSE) |>
    select(subject_id, stim, r, fq) |>
    pivot_wider(names_from = r, values_from = fq) |>
    separate(subject_id, c("subject_id", "sample_id"), sep = "_") |>
    mutate(sample_id = replace_na(sample_id, "1"),
	   sample_id = paste(subject_id, sample_id, sep = "."),
	   stim = recode(stim, "TR7" = "TLR7"),
	   stim = factor(stim, levels = c("unstday0", "BCR", "TLR7", "DN2"))) |>
    arrange(subject_id, sample_id, stim)

write_tsv(meta, 
	  "/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/B_cells_rnaseq/metadata.tsv")


# MGB
batches <- sprintf("04%02d", 1:8)

ids_df <- 
    "/lab-share/IM-Gutierrez-e2/Public/External_datasets/mgb_biobank/%s/etc/ids.tsv" |>
    sprintf(batches) |>
    setNames(batches) |>
    map_df(~read_tsv(., col_names = c("subject_id", "mgb_id")), .id = "batch") |>
    filter(subject_id %in% unique(meta$subject_id))

dir.create("data")
write_tsv(ids_df, "./data/metadata.tsv")

ids_df |>
    select(batch, mgb_id) |>
    group_by(batch) |>
    summarise(id = list(mgb_id)) |>
    ungroup() |>
    mutate(f = paste0("data/", batch, ".txt")) |>
    {function(dat) walk2(dat$id, dat$f, ~write_lines(.x, .y))}()

# Array design for VCF subset
ids_df |>
    distinct(batch) |>
    cross_join(tibble(chr = sprintf("chr%s", 1:22))) |>
    write_tsv("./array_subsetVCF.txt", col_names = FALSE)

# add chrX
ids_df |>
    distinct(batch) |>
    mutate(chr = "chrX") |>
    write_tsv("./array_subsetVCF.txt", col_names = FALSE, append = TRUE)



# for ASE
# need to know the sample id in the VCF
meta |>
    mutate(subject_id = as.numeric(subject_id)) |>
    left_join(ids_df, by = "subject_id") |>
    select(subject_id, sample_id, stim, mgb_id) |>
    write_tsv("./arrayspec_ase.tsv", col_names = FALSE)


