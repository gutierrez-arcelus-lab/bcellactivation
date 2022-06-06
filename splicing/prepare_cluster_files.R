library(tidyverse)

metadata <- 
    "/lab-share/IM-Gutierrez-e2/Public/External_datasets/Scharer/RNAseq/metadata/file_description.tsv" %>%
    read_tsv() %>%
    separate(id, c("status", "dummy", "cellid"), sep = "\\.") %>%
    select(cellid, run, status) %>%
    arrange(cellid, run, status)

# create directories for each cell type
metadata %>%
    pull(cellid) %>%
    unique() %>%
    file.path("./results", .) %>%
    walk(~dir.create(., recursive = TRUE))

# write separate cluster files per cell type
cluster_counts <- "./scharer_perind_numers.counts.gz" %>%
    read.table()

clusters_list <- metadata %>%
    split(.$cellid) %>%
    map(~pull(., run)) %>%
    map(~select(cluster_counts, all_of(.)))


walk(unique(metadata$cellid),
     function(x) {
	 outname = sprintf("./results/%s/%s_perind_numers.counts.gz", x, x)
	 write.table(clusters_list[[x]], file = gzfile(outname), quote = FALSE)
     })

metadata %>%
    split(.$cellid) %>%
    walk(function(x) {
	     cid = unique(x$cellid)
	     outname = sprintf("./results/%s/groups_file.txt", cid)
	     select(x, -cellid) %>% write.table(outname, col.names = FALSE, row.names = FALSE, quote = FALSE)
     })


