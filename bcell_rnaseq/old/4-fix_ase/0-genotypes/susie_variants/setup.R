library(tidyverse)

susie_res <- 
    "../../../../colocalization/finemap/susie_results.tsv" |>
    read_tsv()

susie_fine <- susie_res |>
    filter(!is.na(cs), pip >= 0.1) |>
    distinct(chr, pos)

# Manually include validated variant
val <- 
    tribble(~chr, ~pos,
	    "chr16", 85985027,
	    "chr11", 35073939)

out <- 
    susie_fine |>
    bind_rows(val) |>
    arrange(chr, pos)

write_tsv(out, "susie_variants.tsv", col_names = FALSE)

read_tsv("../array_spec.txt", col_names = c("batch", "chr")) |>
    filter(chr %in% unique(out$chr)) |>
    write_tsv("./array_spec.tsv", col_names = FALSE)
