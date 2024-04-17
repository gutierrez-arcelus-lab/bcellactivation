library(tidyverse)

batches <- sprintf("04%02d", c(1:8, 10)) 
chrs <- paste0("chr", c(1:22, "X"))

# array for subsetting MGB VCFs
array_vcf <- "array_spec_vcf.txt"
array_mbv <- "array_spec_mbv.txt"

expand_grid(batches, chrs) |>
    write_tsv(array_vcf, col_names = FALSE)

# array for running MBV
read_tsv("../0-qc_rnaseq/metadata.tsv") |>
    unite("bam_prefix", c(donor_id, replic_id, stim), sep = "_") |>
    select(bam_prefix) |>
    expand_grid(batch = batches) |>
    write_tsv(array_mbv, col_names = FALSE)
