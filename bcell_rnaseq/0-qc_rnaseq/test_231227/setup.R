library(tidyverse)
library(glue)

labshr <- "/lab-share/IM-Gutierrez-e2/Public"

fastq_files <-
    file.path(labshr, 
	      "Lab_datasets/B_cells_rnaseq", 
	      "resequencing/231227_MG10430_fastq") |>
    list.files(full.names = TRUE, pattern = "\\.fastq\\.gz")

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
    ungroup() |>
    select(-basen, -batch) |>
    pivot_wider(names_from = r, values_from = f) |>
    arrange(donor_id, replic_id, stim)

write_tsv(meta, "metadata.tsv")

# MBV
batches <- sprintf("04%02d", c(1:8, 10))

meta_mbv <-
    meta |>
    unite("prefix", c(donor_id, replic_id, stim), sep = "_") |>
    select(prefix) |>
    expand_grid(batch = batches)

write_tsv(meta_mbv, "./array_spec_mbv.tsv", col_names = FALSE)
