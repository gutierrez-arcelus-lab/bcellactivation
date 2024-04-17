library(tidyverse)

# Create output directory
if (!file.exists("data")) dir.create("data")

mbv <- read_tsv("../../2-mbv/matching_results.tsv")

# Slurm array spec for VCF processing
mbv |>
    distinct(vcf_donor_id, vcf_sample_id, vcf_batch) |>
    cross_join(tibble(chr = sprintf("chr%s", c(1:22, "X")))) |>
    write_tsv("./array_spec.txt", col_names = FALSE)

# meta data to concatenate VCFs for each participant
mbv |>
    distinct(vcf_donor_id, vcf_batch) |>
    write_tsv("./data/vcf_samples.tsv", col_names = FALSE)

