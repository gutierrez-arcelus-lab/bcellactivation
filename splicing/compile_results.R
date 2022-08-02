library(tidyverse)

read_leaf <- function(cell_type, dataset) {
    
    sig <- 
        file.path("./results", dataset, cell_type, "leafcutter_cluster_significance.txt") %>%
        read_tsv() %>%
        separate(cluster, c("chr", "cluster"), sep = ":") %>%
        select(cluster, p, p.adjust, genes)

    eff <- 
        file.path("./results", dataset, cell_type, "leafcutter_effect_sizes.txt") %>%
        read_tsv() %>%
        mutate(cluster = sub("^.*(clu.*)$", "\\1", intron)) %>%
        select(cluster, logef, deltapsi)

    inner_join(sig, eff)
}

cell_types_scharer <- c("rN", "aN", "T3", "SM", "DN") %>%
    setNames(., .)

cell_types_barnas <- "DN"

leaf_scharer_df <- map2_dfr(cell_types_scharer, "scharer", read_leaf, .id = "cell_type") %>%
    group_by(cell_type, cluster) %>%
    slice(which.max(abs(deltapsi))) %>%
    ungroup()

leaf_barnas_df <- map2_dfr(cell_types_barnas, "barnas", read_leaf, .id = "cell_type") %>%
    group_by(cell_type, cluster) %>%
    slice(which.max(abs(deltapsi))) %>%
    ungroup()

leaf_scharer_signif <- leaf_scharer_df %>%
    filter(p.adjust < 0.05, abs(deltapsi) > 0.1, !is.na(genes)) %>%
    separate_rows(genes, sep = ",") %>%
    group_by(genes) %>%
    slice(which.max(abs(deltapsi))) %>%
    ungroup()

leaf_barnas_signif <- leaf_barnas_df %>%
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

write_tsv(leaf_scharer_df, "./results/scharer/leafcutter_filtered.tsv")
write_tsv(leaf_scharer_signif, "./results/scharer/leafcutter_filtered_significant.tsv")
write_tsv(leaf_barnas_df, "./results/barnas/leafcutter_filtered.tsv")
write_tsv(leaf_barnas_signif, "./results/barnas/leafcutter_filtered_significant.tsv")


write_tsv(sqtl_df, "./mu2021_data/filtered_sqtl.tsv")
