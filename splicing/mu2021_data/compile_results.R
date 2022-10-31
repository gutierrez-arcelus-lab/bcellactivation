library(tidyverse)


sqtl_df <- list.files(".", pattern = "txt\\.gz") %>%
    setNames(., sub("\\.txt\\.gz", "", basename(.))) %>%
    map_df(. %>% 
	   read_tsv(col_types = c(sid = "c")) %>%
	   select(pid, genes, sid, slope, qval) %>%
	   filter(qval < 0.05, !is.na(genes)) %>%
	   separate_rows(genes, sep = ",") %>%
	   group_by(genes) %>%
	   slice(which.min(qval)) %>%
	   ungroup(),
       .id = "dataset")

write_tsv(sqtl_df, "./filtered_sqtl.tsv")
