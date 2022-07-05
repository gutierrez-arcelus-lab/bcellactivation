library(tidyverse)
library(leafcutter)

exon_table <- read_tsv("./hg19_exon.txt.gz")

sqtl_dir <- "/lab-share/IM-Gutierrez-e2/Public/External_datasets/Mu2021/eqtl_sqtl_summ_stats"

qtl_df <- list.files(sqtl_dir, recursive = TRUE, full.names = TRUE) %>%
    grep("sQTL", ., value = TRUE) %>%
    setNames(., str_remove(., sqtl_dir)) %>%
    setNames(., str_remove(names(.), "^/")) %>%
    setNames(., str_remove(names(.), "\\.txt\\.gz$")) %>%
    map_df(~read_tsv(., col_types = "cidddcdddddd") %>%
	   mutate(pid = str_replace(pid, ":(?=[-+])", "_"),
		  pid = paste0("chr", pid)), 
	   .id = "source") %>%
    separate(source, c("source", "dataset"), sep = "/") %>%
    mutate(source = str_remove(source, "_sQTL"),
	   dataset = ifelse(is.na(dataset) & source == "DGN", "Whole Blood", dataset))

intron_meta <- get_intron_meta(qtl_df$pid)

clu_gene_map <- map_clusters_to_genes(intron_meta, exon_table) %>%
    as_tibble()

qtl_annot <- qtl_df %>%
    mutate(clu = sub("^(chr[^:]+):[^:]+:[^:]+:(.+)$", "\\1:\\2", pid)) %>%
    left_join(clu_gene_map, by = "clu") %>%
    select(-clu)

write_tsv(qtl_annot, "./mu2021_data/sqtl_all.tsv")
