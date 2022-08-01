library(tidyverse)

read_leaf <- function(cell_type) {
    
    sig <- 
        file.path("./results", cell_type, "leafcutter_cluster_significance.txt") %>%
        read_tsv() %>%
        separate(cluster, c("chr", "cluster"), sep = ":") %>%
        select(cluster, p, p.adjust, genes)

    eff <- 
        file.path("./results", cell_type, "leafcutter_effect_sizes.txt") %>%
        read_tsv() %>%
        mutate(cluster = sub("^.*(clu.*)$", "\\1", intron)) %>%
        select(cluster, logef, deltapsi)

    inner_join(sig, eff)
}

gene_strand <- 
    "/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens/gencode.v38.primary_assembly.annotation.gtf" %>%
    read_tsv(comment = "#", col_names = FALSE) %>%
    filter(X3 == "gene") %>%
    mutate(gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+")) %>%
    select(gene_name, strand = X7)

cell_types <- c("rN", "aN", "T3", "SM", "DN") %>%
    setNames(., .)

leaf_df <- map_df(cell_types, read_leaf, .id = "cell_type") %>%
    group_by(cell_type, cluster) %>%
    slice(which.max(abs(deltapsi))) %>%
    ungroup()


# filter according to gene strand orientation
# recompute adjusted p-values
leaf_filt_df <- leaf_df %>%
  mutate(strand_clu = str_extract(cluster, "([+-]$)")) %>%
  separate_rows(genes, sep = ",") %>%
  left_join(gene_strand, by = c("genes" = "gene_name")) %>%
  group_by(cell_type, cluster) %>%
  filter(n_distinct(strand) == 1) %>%
  ungroup() %>%
  filter(strand_clu == strand) %>%
  group_by(cell_type, cluster, p, p.adjust, logef, deltapsi) %>%
  summarise(genes = paste(genes, collapse = ",")) %>%
  group_by(cell_type) %>%
  mutate(p.adjust = p.adjust(p, "fdr")) %>%
  ungroup()

leaf_signif <- leaf_filt_df %>%
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

write_tsv(leaf_filt_df, "./results/leafcutter_filtered.tsv")
write_tsv(leaf_signif, "./results/leafcutter_filtered_significant.tsv")
write_tsv(sqtl_df, "./mu2021_data/filtered_sqtl.tsv")
