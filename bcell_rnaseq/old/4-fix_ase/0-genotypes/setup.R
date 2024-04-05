library(tidyverse)

# Create output directory
if (!file.exists("data")) dir.create("data")

mbv <- read_tsv("../../3-mbv/match_summary.tsv")

mbv |>
    distinct(vcf_batch, vcf_id) |>
    group_by(vcf_batch) |>
    summarise(id = list(vcf_id)) |>
    ungroup() |>
    mutate(f = paste0("data/", vcf_batch, ".txt")) |>
    {function(dat) walk2(dat$id, dat$f, ~write_lines(.x, .y))}()

# Slurm array spec for VCF processing
mbv |>
    distinct(vcf_batch) |>
    cross_join(tibble(chr = sprintf("chr%s", c(1:22, "X")))) |>
    write_tsv("./array_spec.txt", col_names = FALSE)

# Metadata for creating individual VCFs
mbv |>
    distinct(vcf_donor_id, vcf_id) |>
    write_tsv("./data/metadata.tsv", col_names = FALSE)
