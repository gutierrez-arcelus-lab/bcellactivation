library(tidyverse)

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

batches <- sprintf("04%02d", c(1:8, 10))
chrs <- c(1:22, "X")

var_df <- expand_grid(batch = batches, chr = chrs) %>%
    mutate(path = map2_chr(chr, batch, ~sprintf("/temp_work/ch229163/chr%s_%s_typed.txt", .x, .y)),
	   data = map(path, read_lines))

out_df <- var_df %>%
    select(-path) %>%
    unnest(c(data)) %>%
    group_split(batch) %>%
    map(~select(., -batch)) %>%
    reduce(inner_join, by = c("chr", "data")) %>%
    group_nest(chr) %>%
    mutate(data = map(data, unlist),
	   out = sprintf("/temp_work/ch229163/chr%s_typed.txt", chr))

out_df %>%
    walk2(.x = .$data, .y = .$out, .f = ~write_lines(.x, .y))
