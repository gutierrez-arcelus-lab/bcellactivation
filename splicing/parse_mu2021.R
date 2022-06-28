library(tidyverse)
library(leafcutter)

exon_table <- read_tsv("./hg19_exon.txt.gz")

dataset <- commandArgs(TRUE)[1]

qtl <- 
    file.path("/lab-share/IM-Gutierrez-e2/Public/External_datasets/Mu2021/eqtl_sqtl_summ_stats", dataset) %>%
    read_tsv(col_types = "cidddcdddddd") %>%
    mutate(pid = str_replace(pid, ":(?=[-+])", "_"),
	   pid = paste0("chr", pid))

intron_meta <- get_intron_meta(qtl$pid)

clu_gene_map <- map_clusters_to_genes(intron_meta, exon_table) %>%
    as_tibble()

qtl_annot <- qtl %>%
    mutate(clu = sub("^(chr[^:]+):[^:]+:[^:]+:(.+)$", "\\1:\\2", pid)) %>%
    left_join(clu_gene_map, by = "clu") %>%
    select(-clu)

out <- dataset %>%
    tolower() %>%
    sub("/", "_", .)

write_tsv(qtl_annot, file.path("./mu2021_data", out))
