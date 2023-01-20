library(tidyverse)

fqdir <- "/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/B_cells_rnaseq/34.198.31.178/221025_MG10430_fastq/"

meta <- tibble(fq = list.files(fqdir)) |>
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
