library(tidyverse)

meta <- read_tsv("./mbv.spec", col_names = FALSE)

files <- file.path("results", paste0(meta$X1, "_", meta$X2, ".txt"))
names(files) <- meta$X1

res <- files |> 
    map_dfr(~read_delim(., delim = " "), .id = "bamid") |>
    select(bamid, vcfid = SampleID, perc_het_consistent, perc_hom_consistent)

matches_df <- res |> 
    group_by(bamid) |>
    filter(perc_het_consistent == max(perc_het_consistent)) |>
    ungroup()

res |> 
    filter(! vcfid %in% unique(matches_df$vcfid) ) |>
    group_by(bamid) |>
    filter(perc_het_consistent == max(perc_het_consistent)) |>
    ungroup()



