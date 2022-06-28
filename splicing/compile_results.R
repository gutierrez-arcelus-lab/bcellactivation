library(tidyverse)

read_leaf <- function(cell_type) {
    
    sig <- 
        file.path("./results", cell_type, "leafcutter_cluster_significance.txt") %>%
        read_tsv() %>%
        separate(cluster, c("chr", "cluster"), sep = ":") %>%
        select(cluster, p.adjust, genes)

    eff <- 
        file.path("./results", cell_type, "leafcutter_effect_sizes.txt") %>%
        read_tsv() %>%
        mutate(cluster = sub("^.*(clu.*)$", "\\1", intron)) %>%
        select(cluster, logef, deltapsi)

    inner_join(sig, eff)
}

cell_types <- c("rN", "aN", "T3", "SM", "DN") %>%
    setNames(., .)

leaf_df <- map_df(cell_types, read_leaf, .id = "cell_type") %>%
    group_by(cell_type, cluster) %>%
    slice(which.max(abs(deltapsi))) %>%
    ungroup()

leaf_filtered <- leaf_df %>%
    filter(p.adjust < 0.05, abs(deltapsi) > 0.1, !is.na(genes)) %>%
    separate_rows(genes, sep = ",") %>%
    group_by(genes) %>%
    slice(which.max(abs(deltapsi))) %>%
    ungroup()

sqtl_df <- list.files("./mu2021_data", full.names = TRUE) %>%
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

write_tsv(leaf_filtered, "./results/filtered_leafcutter.tsv")
write_tsv(sqtl_df, "./mu2021_data/filtered_sqtl.tsv")

