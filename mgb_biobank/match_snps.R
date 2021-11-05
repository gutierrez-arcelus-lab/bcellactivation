library(tidyverse)

vars_path <- commandArgs(TRUE)[1]
vars_1kg_path <- commandArgs(TRUE)[2]
out_path <- commandArgs(TRUE)[3] 

read_variants <- function(var_path) {

    var_path %>%
	read_delim(delim = " ") %>%
	setNames(sub("#", "", names(.))) %>%
	add_count(POS) %>%
	filter(n == 1) %>%
	filter(!(REF == "A" & ALT == "T"),
	       !(REF == "T" & ALT == "A"),
	       !(REF == "C" & ALT == "G"),
	       !(REF == "G" & ALT == "C"))
}

vars_df <- read_variants(vars_path)

vars_1000G <- read_variants(vars_1kg_path) %>%
    mutate(CHROM = ifelse(is.numeric(CHROM), paste0("chr", CHROM), CHROM))

merge_df <- inner_join(vars_df, vars_1000G, by = c("CHROM", "POS", "REF", "ALT"))

out <- merge_df %>%
    select(1:2) %>%
    arrange(POS)

write_tsv(out, out_path, col_names = FALSE)
