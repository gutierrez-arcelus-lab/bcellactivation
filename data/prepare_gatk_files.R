library(tidyverse)

report <- "./GCF_000001405.39_GRCh38.p13_assembly_report.txt" %>%
    read_tsv(comment = "#", col_names = FALSE) %>%
    select(X1, X2, X7, X8, X10)

report %>%
    select(X7, X10) %>%
    mutate(X10 = sub("^chr[^_]+_", "", X10),
	   X10 = sub("_random$", "", X10),
	   X10 = case_when(grepl("v[12]$", X10) ~ sub("v", ".", X10),
			   TRUE ~ X10)) %>%
    write_tsv("chr_names.txt", col_names = FALSE)

report %>%
    filter(X8 == "Primary Assembly") %>%
    filter(X7 != "na") %>%
    pull(X7) %>%
    write_lines("pri_chrs.txt")
