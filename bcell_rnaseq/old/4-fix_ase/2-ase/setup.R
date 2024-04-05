library(tidyverse)

meta <- 
    file.path("/lab-share/IM-Gutierrez-e2/Public",
	      "Lab_datasets/B_cells_rnaseq/metadata.tsv") |>
    read_tsv()

mbv <- read_tsv("../../3-mbv/match_summary.tsv")

fixed_meta <- 
    left_join(meta, mbv, join_by(donor_id == bam_donor_id, sample_id, stim == bam_stim)) |>
    mutate(rep_id = sub("^[^\\.]+\\.(\\d)$", "\\1", sample_id)) |>
    mutate(sample_id = paste(vcf_donor_id, rep_id, sep = ".")) |>
    select(donor_id = vcf_donor_id, sample_id, stim, vcf_id)

write_tsv(fixed_meta, "./array_spec.tsv", col_names = FALSE)
