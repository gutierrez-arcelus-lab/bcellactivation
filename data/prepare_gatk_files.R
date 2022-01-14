library(tidyverse)

report <- commandArgs(TRUE)[1] %>%
    read_tsv(comment = "#", col_names = FALSE) %>%
    filter(X7 != "na", X10 != "na") %>%
    select(X1, X2, X7, X8, X10) %>%
    mutate(genc = sub("^chr[^_]+_", "", X10),
	   genc = sub("_.+$", "", genc),
	   genc = case_when(grepl("v[12]$", genc) ~ sub("v", ".", genc),
			   TRUE ~ genc))
report %>%
    select(X7, X10) %>%
    write_tsv("/temp_work/ch229163/chr_dbsnpToGencode_names.txt", col_names = FALSE)

report %>%
    select(X10, genc) %>%
    write_tsv("/temp_work/ch229163/chr_UcscToGencode_names.txt", col_names = FALSE)

report %>%
    filter(X8 == "Primary Assembly") %>%
    select(genc) %>%
    mutate(start = 0, end = 1e9) %>%
    write_tsv("ref_pri_chr.bed", col_names = FALSE)
