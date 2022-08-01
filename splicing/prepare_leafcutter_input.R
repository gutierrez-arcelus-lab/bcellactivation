library(tidyverse)

# Scharer
metadata <- 
    "/lab-share/IM-Gutierrez-e2/Public/External_datasets/Scharer/RNAseq/metadata/file_description.tsv" %>%
    read_tsv() %>%
    separate(id, c("status", "dummy", "cellid"), sep = "\\.") %>%
    select(cellid, run, status) %>%
    arrange(cellid, run, status)

## create directories for each cell type
metadata %>%
    pull(cellid) %>%
    unique() %>%
    file.path("./results/scharer", .) %>%
    walk(~dir.create(., recursive = TRUE))

## junction files
juncfiles <- 
    "/lab-share/IM-Gutierrez-e2/Public/vitor/lupus/splicing/mapping/scharer" %>%
    list.files(pattern = "*.junc$", full.names = TRUE)

tibble(f = juncfiles, run = str_extract(juncfiles, "SRR\\d+")) %>%
    left_join(metadata) %>%
    group_by(cellid) %>%
    nest() %>%
    ungroup() %>%
    {walk2(.$cellid, .$data, 
           ~write_lines(.y$f, file.path("./results/scharer", .x, "juncfiles.txt")))}

## groups file
metadata %>%
    group_by(cellid) %>%
    nest() %>%
    {walk2(.$cellid, .$data, 
           ~write_tsv(.y, file.path("./results/scharer", .x, "groups_file.txt"), col_names = FALSE))}

# Barnas
metadata_barnas <- 
    "/lab-share/IM-Gutierrez-e2/Public/External_datasets/Barnas/SraRunTable.txt" %>%
    read_csv() %>%
    mutate(cellid = case_when(flow_sort_selection == "CD19+ CD3- CD20+ IgD- CD27-" ~ "DN",
			      flow_sort_selection == "CD19+ CD3- CD20+ IgD+ CD27-" ~ "Naive B",
			      flow_sort_selection == "CD19- CD14+ CD3-" ~ "Monocytes")) %>%
    filter(cellid == "DN") %>%
    select(cellid, run = Run, status = disease_state) 

dir.create("./results/barnas/DN")

##
juncfiles_barnas <- 
    "/lab-share/IM-Gutierrez-e2/Public/vitor/lupus/splicing/mapping/barnas" %>%
    list.files(pattern = "*.junc$", full.names = TRUE)

tibble(f = juncfiles_barnas, run = str_extract(juncfiles_barnas, "SRR\\d+")) %>%
    left_join(metadata_barnas) %>%
    group_by(cellid) %>%
    nest() %>%
    ungroup() %>%
    {walk2(.$cellid, .$data, 
           ~write_lines(.y$f, file.path("./results/barnas", .x, "juncfiles.txt")))}

## groups file
metadata_barnas %>%
    arrange(cellid, status, run) %>%
    group_by(cellid) %>%
    nest() %>%
    {walk2(.$cellid, .$data, 
           ~write_tsv(.y, file.path("./results/barnas", .x, "groups_file.txt"), col_names = FALSE))}
    
