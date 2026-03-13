library(tidyverse)
library(glue)

stims <- read_lines("./stims.txt")

results <-
    glue("./results_vert/{stims}/knownResults.txt") |>
    #glue("../atacseq/homer/results_vert_mRp/{stims}/knownResults.txt") |>
    setNames(paste0(stims, "c")) |>
    map_dfr(~read_tsv(.) |>
	    janitor::clean_names() |>
	    {function(x) setNames(x, str_remove(names(x), "_of_\\d+$"))}(),
	    .id = "stim")

write_tsv(results, "./results_vert/results.tsv")
