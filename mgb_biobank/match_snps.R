library(tidyverse)

CHR <- commandArgs(TRUE)[1]

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

vars_df <- sprintf("/temp_work/ch229163/chr%s.variants.txt.gz", CHR) %>%
    read_variants()

vars_1000G <- sprintf("/temp_work/ch229163/chr%s.1000G.variants.txt.gz", CHR) %>%
    read_variants()

merge_df <- inner_join(vars_df, vars_1000G, by = c("CHROM", "POS", "REF", "ALT"))

out <- merge_df %>%
    select(1:2) %>%
    arrange(POS)

write_tsv(out, 
	  sprintf("/temp_work/ch229163/chr%s.filtered.pos", CHR),
	  col_names = FALSE)
