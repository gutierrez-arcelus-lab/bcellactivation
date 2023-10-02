library(tidyverse)

samples <- read_lines("./samples.txt") %>%
    c(read_lines("./samples_hc.txt"))

quants_df <- file.path("./salmon", samples, "quant.sf") %>%
    setNames(samples) %>%
    map_df(read_tsv, .id = "sample_id")

quants_wide <- quants_df %>%
    filter(!grepl("^HHV4", Name)) %>%
    select(sample_id, transcript_id = Name, TPM) %>%
    pivot_wider(names_from = sample_id, values_from = TPM) 

write_tsv(quants_wide, "./compiled_expression_humangenes.tsv")
