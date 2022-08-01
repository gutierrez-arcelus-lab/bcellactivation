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

# junction files
juncfiles <- 
    "/lab-share/IM-Gutierrez-e2/Public/vitor/lupus/splicing/mapping" %>%
    list.files(pattern = "*.junc$", full.names = TRUE)

tibble(f = juncfiles, run = str_extract(juncfiles, "SRR\\d+")) %>%
    left_join(metadata) %>%
    group_by(cellid) %>%
    nest() %>%
    ungroup() %>%
    {walk2(.$cellid, .$data, 
           ~write_lines(.y$f, file.path("./results", .x, "juncfiles.txt")))}

# groups file
metadata %>%
    group_by(cellid) %>%
    nest() %>%
    {walk2(.$cellid, .$data, 
           ~write_tsv(.y, file.path("./results", .x, "groups_file.txt"), col_names = FALSE))}
