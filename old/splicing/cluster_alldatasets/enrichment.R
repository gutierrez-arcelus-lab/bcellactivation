library(fgsea)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)

run_enrichment <- function(gene_list) {

    enrichGO(gene = gene_list,
	OrgDb = org.Hs.eg.db,
	ont = "BP",
	keyType = "ENSEMBL",
	pAdjustMethod = "fdr",
	pvalueCutoff = 0.1,
	qvalueCutoff = 0.1,
	readable = TRUE) %>%
    as.data.frame() %>%
    as_tibble()
}

# Data
gene_id_df <- 
    "/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens/gencode.v38.txToGene.tsv" %>%
    read_tsv() %>%
    distinct(gene_id, gene_name) %>%
    filter(!grepl("PAR_Y$", gene_id)) %>%
    mutate(gene_id = str_remove(gene_id, "\\.\\d+$")) %>%
    add_count(gene_name) %>%
    filter(n == 1) %>%
    select(-n)

cell_types <- c("aN", "DN", "rN", "T3", "SM")

leaf_df <- 
    read_tsv("results/scharer/leafcutter_filtered.tsv") %>%
    filter(!is.na(genes)) %>%
    separate_rows(genes, sep = ",") %>%
    group_by(cell_type, genes) %>%
    slice(which.max(absdpsi)) %>%
    ungroup() %>%
    inner_join(gene_id_df, by = c("genes" = "gene_name"))

go_res <- leaf_df %>%
    filter(p.adjust < .1) %>%
    split(.$cell_type) %>%
    map("gene_id") %>%
    map_df(run_enrichment, .id = "cell_type")

write_tsv(go_res, "./results/scharer/go.txt")

# GSEA
hallmark <- gmtPathways("/lab-share/IM-Gutierrez-e2/Public/References/msigdb/h.all.v2022.1.Hs.symbols.gmt.txt")
pathwaysgo <- gmtPathways("/lab-share/IM-Gutierrez-e2/Public/References/msigdb/c5.all.v2022.1.Hs.symbols.gmt.txt")
gobp <- keep(pathwaysgo, grepl("^GOBP", names(pathwaysgo)))
humanphenos <- keep(pathwaysgo, grepl("^HP", names(pathwaysgo)))

all_genes <- unlist(c(pathwaysgo, hallmark))

gsea_list <- leaf_df %>%
    filter(genes %in% all_genes) %>%
    select(cell_type, SYMBOL = genes, stat = absdpsi) %>%
    arrange(cell_type, stat) %>%
    split(.$cell_type) %>%
    map(~select(., -cell_type)) %>%
    map(deframe)

scharer_gsea_hallmark_res <- 
    map(gsea_list, ~fgsea(pathways = hallmark, stats = ., scoreType = "pos"))

bind_rows(scharer_gsea_hallmark_res, .id = "cell_type") %>%
    as_tibble() %>%
    filter(padj < 0.2)

scharer_gsea_gobp_res <- 
    map(gsea_list, ~fgsea(pathways = gobp, stats = ., scoreType = "pos"))

bind_rows(scharer_gsea_gobp_res, .id = "cell_type") %>%
    as_tibble() %>%
    group_by(cell_type) %>%
    slice_max(n = 5, abs(NES)) %>%
    ungroup() %>% print(n = Inf)

scharer_gsea_phenos_res <- 
    map(gsea_list, ~fgsea(pathways = humanphenos, stats = ., scoreType = "pos"))

bind_rows(scharer_gsea_phenos_res, .id = "cell_type") %>%
    as_tibble() %>%
    arrange(pval)
