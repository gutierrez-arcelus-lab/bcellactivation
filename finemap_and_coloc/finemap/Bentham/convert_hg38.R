library(tidyverse)

# Try to convert using the GRCh38 coordinates for the same variants reported
# in the summary statistics harmonized by GWAS catalog

gwas_hg38 <- 
    "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Bentham/harmonized/26502338-GCST003156-EFO_0002690.h.tsv.gz" |>
    data.table::fread() |>
    as_tibble() |>
    select(rsid = variant_id, chrom = chromosome, pos = base_pair_location, other_allele, effect_allele)

gwas_hg37 <- read_tsv("./data/summary_stats.tsv")

out <- 
    left_join(gwas_hg37, gwas_hg38, join_by(chrom, rsid)) |>
    filter(!is.na(pos.y)) |>
    distinct(locus, chrom, pos = pos.y, rsid, ref, alt, beta, se, logp)

write_tsv(out, "./data/summary_stats_hg38.tsv")
