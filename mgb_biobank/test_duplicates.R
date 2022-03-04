library(tidyverse)

ids_df <- "/lab-share/IM-Gutierrez-e2/Public/mgb_biobank/%04d/etc/ids.tsv" %>%
    sprintf(c(401:408, 410)) %>%
    setNames(c(paste0("0", c(401:408, 410)))) %>%
    map_df(read_tsv, col_names = FALSE, .id = "batch") %>%
    select(subject_id = X1, sample_id = X2, batch)

dups_df <- ids_df %>%
    add_count(subject_id) %>%
    filter(n > 1) %>%
    arrange(subject_id, sample_id)

dups_df %>%
    select(subject_id, batch) %>%
    arrange(subject_id, batch) %>%
    group_by(subject_id) %>%
    mutate(i = 1:n()) %>%
    ungroup() %>%
    pivot_wider(names_from = i, values_from = batch) %>%
    count(`1`, `2`)


batches <- sprintf("04%02d", c(1:8, 10))

females <- "./results/females_%s.txt" %>%
    sprintf(batches) %>%
    setNames(batches) %>%
    map_df(~tibble(sampleid = read_lines(.)), .id = "batch")

dups_df %>%
    mutate(female = ifelse(sample_id %in% females$sampleid, 1, 0)) %>%
    distinct(subject_id, .keep_all = TRUE) %>%
    summarise(mean(female))

ids_df %>%
    filter(! subject_id %in% dups_df$subject_id) %>%
    mutate(female = ifelse(sample_id %in% females$sampleid, 1, 0)) %>%
    summarise(mean(female))

dups_df %>%
    select(sample_id, batch) %>%
    arrange(batch, sample_id) %>%
    group_nest(batch) %>%
    mutate(out = paste0("./results/duplicates_", batch, ".txt"),
	   data = map(data, unlist)) %>%
    walk2(.x = .$data, .y = .$out, .f = ~write_lines(.x, .y))
