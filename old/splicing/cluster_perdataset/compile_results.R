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

bind_rows(scharer = leaf_scharer_df, barnas = leaf_barnas_df, andreoletti = leaf_and_df, .id = "dataset") %>%
    select(dataset, cell_type, genes, p, p.adjust, absdpsi) %>%
    write_tsv("./results/leafcutter_merged_filtered.tsv")

bind_rows(scharer = leaf_scharer_df, andreoletti = leaf_and_df, .id = "dataset") %>%
    select(dataset, cell_type, genes, cluster, p, p.adjust, absdpsi) %>%
    mutate_at(vars(p, p.adjust), ~format(., scientific = FALSE)) %>%
    write_tsv("./results/leafcutter_ScharerAndreoletti_filt.tsv")

####
cell_type <- "DN"
dataset <- "scharer"

read_leaf_raw <- function(cell_type, dataset) {
    
    sig <- 
        file.path("./results", dataset, cell_type, "leafcutter_cluster_significance.txt") %>%
        read_tsv() %>%
        separate(cluster, c("chr", "cluster"), sep = ":") %>%
        select(cluster, p, p.adjust, genes) %>%
	filter(!is.na(p))

    eff <- 
        file.path("./results", dataset, cell_type, "leafcutter_effect_sizes.txt") %>%
        read_tsv() %>%
	extract(intron, c("junction", "cluster"), "([^:]+:\\d+:\\d+):(clu_.+)") %>%
        select(cluster, junction, deltapsi)

    inner_join(sig, eff)
}

scharer_all <- map2_dfr(cell_types_scharer, "scharer", read_leaf_raw, .id = "cell_type")

write_tsv(scharer_all, "./results/scharer/leafcutter_all_significant.tsv")
