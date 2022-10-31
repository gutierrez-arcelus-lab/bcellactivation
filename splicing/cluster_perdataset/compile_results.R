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

cell_types_barnas <- c("DN", "Naive") %>%
    setNames(., .)

cell_types_and <- setNames("B", "B")

leaf_scharer_df <- map2_dfr(cell_types_scharer, "scharer", read_leaf, .id = "cell_type") %>%
    mutate(absdpsi = abs(deltapsi)) %>%
    group_by(cell_type, cluster) %>%
    slice(which.max(absdpsi)) %>%
    ungroup() %>%
    select(cell_type, cluster, p, p.adjust, genes, logef, absdpsi) %>%
    arrange(cell_type, p)

leaf_barnas_df <- map2_dfr(cell_types_barnas, "barnas", read_leaf, .id = "cell_type") %>%
    mutate(absdpsi = abs(deltapsi)) %>%
    group_by(cell_type, cluster) %>%
    slice(which.max(absdpsi)) %>%
    ungroup() %>%
    select(cell_type, cluster, p, p.adjust, genes, logef, absdpsi) %>%
    arrange(cell_type, p)

leaf_and_df <- map2_dfr(cell_types_and, "andreoletti", read_leaf, .id = "cell_type") %>%
    mutate(absdpsi = abs(deltapsi)) %>%
    group_by(cell_type, cluster) %>%
    slice(which.max(absdpsi)) %>%
    ungroup() %>%
    select(cell_type, cluster, p, p.adjust, genes, logef, absdpsi) %>%
    arrange(cell_type, p)

write_tsv(leaf_scharer_df, "./results/scharer/leafcutter_filtered.tsv")
write_tsv(leaf_barnas_df, "./results/barnas/leafcutter_filtered.tsv")
write_tsv(leaf_and_df, "./results/andreoletti/leafcutter_filtered.tsv")
