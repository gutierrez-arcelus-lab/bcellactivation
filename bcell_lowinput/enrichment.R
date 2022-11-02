library(fgsea)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)

cpm_df <- read_tsv("./data/edger_cpm_fit.tsv")
edger_results <- read_tsv("./data/edger_de_genes.tsv")

pathwaysgo <- 
    "/lab-share/IM-Gutierrez-e2/Public/References/msigdb/c5.all.v2022.1.Hs.symbols.gmt.txt" %>%
    gmtPathways()

humanphenos <- keep(pathwaysgo, grepl("^HP", names(pathwaysgo)))

edger_filt <- edger_results %>%
    group_by(stim, gene_name) %>%
    add_count() %>%
    ungroup() %>%
    filter(n == 1) 

gene_list <- edger_filt %>%
    select(stim, SYMBOL = gene_name, stat = PValue) %>%
    arrange(stim, stat) %>%
    split(.$stim) %>%
    map(~select(., -stim)) %>%
    map(deframe)

gsea_res <- gene_list %>%
    map_df(~fgsea(pathways = humanphenos, stat = ., scoreType = "pos"), .id = "stim") %>%
    as_tibble()

write_tsv(gsea_res, "./data/gsea_results.tsv")
