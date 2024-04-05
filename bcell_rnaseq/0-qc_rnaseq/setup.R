library(tidyverse)

labshr <- "/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/B_cells_rnaseq"

fastq_files <-
    c("34.198.31.178/221025_MG10430_fastq", 
      "resequencing/231227_MG10430_fastq",
      "resequencing/240122_MG10430_fastq",
      "resequencing/240320_MG10430_fastq") |>
    map_chr(~file.path(labshr, .)) |>
    map(~list.files(., full.names = TRUE, pattern = "\\.fastq\\.gz")) |>
    unlist()

meta <- 
    tibble(f = fastq_files) |>
    mutate(basen = basename(f),
	   r = str_extract(basen, "_(R[12])_", group = 1)) |>
    extract(basen, 
	    c("batch", "dummy", "donor_id", "replic_id", "stim"),
	    "(\\d+)_(MG_)?(\\d+)[_-](\\d[_-])?([^_-]+)",
	    remove = FALSE) |>
    select(batch, donor_id, replic_id, stim, f, basen, r) |>
    mutate(replic_id = parse_number(replic_id),
	   replic_id = replace_na(replic_id, 1),
	   stim = recode(stim, "TR7" = "TLR7")) |>
    group_by(batch, donor_id, replic_id, stim) |>
    filter(all(c("R1", "R2") %in% r)) |>
    group_by(donor_id, replic_id, stim, r) |>
    summarise(fastq = paste(f, collapse = ",")) |>
    ungroup() |>
    pivot_wider(names_from = r, values_from = fastq)

write_tsv(meta, "metadata.tsv")
