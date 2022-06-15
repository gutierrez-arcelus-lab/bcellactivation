library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)

select <- dplyr::select
slice <- dplyr::slice

gene_id_df <- 
    "/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens/gencode.v38.txToGene.tsv" %>%
    read_tsv() %>%
    distinct(gene_id, gene_name) %>%
    filter(!grepl("PAR_Y$", gene_id)) %>%
    mutate(gene_id = str_remove(gene_id, "\\.\\d+$"))

# Leafcutter results
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

run_enrichment <- function(gene_list) {

    enrichGO(gene = gene_list,
	OrgDb = org.Hs.eg.db,
	ont = "BP",
	keyType = "ENSEMBL",
	pAdjustMethod = "fdr",
	pvalueCutoff = 0.05,
	qvalueCutoff = 0.05,
	readable = TRUE) %>%
    as.data.frame() %>%
    as_tibble()
}

cell_types <- list.files("results") %>%
    .[!. %in% "PbPc"] %>%
    setNames(., .)

leaf_df <- map_df(cell_types, read_leaf, .id = "cell_type") %>%
    group_by(cell_type, cluster) %>%
    slice(which.min(p.adjust)) %>%
    ungroup() %>%
    filter(p.adjust < 0.05, abs(deltapsi) >= 0.1)

leaf_genes <- leaf_df %>%
    select(cell_type, genes) %>%
    filter(!is.na(genes)) %>%
    separate_rows(genes, sep = ",") %>%
    distinct(cell_type, genes) %>%
    left_join(gene_id_df, by = c("genes" = "gene_name"))

go_res <- leaf_genes %>%
    split(.$cell_type) %>%
    map("gene_id") %>%
    map_df(run_enrichment, .id = "cell_type")

write_tsv(go_res, "./results/enrichment.tsv")
