library(tidyverse)

# Gene Ontology analysis
library(clusterProfiler)
library(org.Hs.eg.db)

count <- dplyr::count
select <- dplyr::select
slice <- dplyr::slice

deg <- read_tsv("../bcell_lowinput/results/edger/diff_expr_all_times.tsv")

# GO
run_enrichment <- function(gene_list) {

    enrichGO(gene = gene_list,
	OrgDb = org.Hs.eg.db,
	ont = "BP",
	keyType = "ENSEMBL",
	pAdjustMethod = "fdr",
	qvalueCutoff = 0.1,
	readable = TRUE)
}

deg_genes <- 
    deg |> 
    filter(group1 == "BCR.72", group2 == "DN2.72") |>
    filter(logFC > 0, FDR <= 0.01) |>
    mutate(gene_id = str_remove(gene_id, "\\.\\d+$")) |>
    pull(gene_id)

go_res <- 
    run_enrichment(deg_genes) |>
    as_tibble()

write_tsv(go_res, "./table-s1.tsv")
