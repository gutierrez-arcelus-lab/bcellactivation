library(tidyverse)
library(glue)

stims <- read_lines("../03_atacseq/3-motif_analysis/data/stims.txt")

results <-
    glue("../03_atacseq/3-motif_analysis/results/{stims}/knownResults.txt") |>
    setNames(paste0(stims, "c")) |>
    map_dfr(~read_tsv(.) |>
            janitor::clean_names() |>
            {function(x) setNames(x, str_remove(names(x), "_of_\\d+$"))}(),
            .id = "stim")

write_tsv(results, "./data/supplementary_data_3_homer.tsv")
