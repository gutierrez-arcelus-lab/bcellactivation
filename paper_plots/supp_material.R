library(tidyverse)

# Gene Ontology analysis
library(clusterProfiler)
library(org.Hs.eg.db)

count <- dplyr::count
select <- dplyr::select
slice <- dplyr::slice

deg <- read_tsv("../bcell_lowinput/results/edger/diff_expr_all_times_all_genes.tsv")

contrasts_df <-
    deg |>
    distinct(contrast) |>
    separate(contrast, c("stim_test", "stim_base"), sep = "-", remove = FALSE)


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
    left_join(contrasts_df, join_by(contrast)) |>
    filter(stim_base == "BCR.72", stim_test == "DN2.72") |>
    filter(logFC > 0, FDR <= 0.01) |>
    mutate(gene_id = str_remove(gene_id, "\\.\\d+$")) |>
    pull(gene_id)

go_res <- 
    run_enrichment(deg_genes) |>
    as_tibble()

write_tsv(go_res, "./Supplementary_Table_4.tsv")
