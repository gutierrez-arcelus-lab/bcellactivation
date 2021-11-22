library(tidyverse)

paths <- commandArgs(TRUE)

var_paths <- paths[-length(paths)]
out_path <- last(paths)

read_variants <- function(var_path) {

    var_path %>%
	read_delim(delim = " ") %>%
	setNames(sub("#", "", names(.))) %>%
	add_count(POS) %>%
	filter(n == 1) %>%
	filter(!(REF == "A" & ALT == "T"),
	       !(REF == "T" & ALT == "A"),
	       !(REF == "C" & ALT == "G"),
	       !(REF == "G" & ALT == "C")) %>%
	select(-n)
}

out_df <- c(mgb_paths, kgp_path) %>%
    map(read_variants) %>%
    reduce(inner_join, by = c("CHROM", "POS", "REF", "ALT")) %>%
    select(CHROM, POS) %>%
    arrange(POS)

write_tsv(out_df, out_path, col_names = FALSE)
