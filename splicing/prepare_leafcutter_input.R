library(tidyverse)

# Scharer
metadata_scharer <- 
    "/lab-share/IM-Gutierrez-e2/Public/External_datasets/Scharer/RNAseq/metadata/file_description.tsv" %>%
    read_tsv() %>%
    separate(id, c("status", "dummy", "cellid"), sep = "\\.") %>%
    select(cellid, run, status) %>%
    arrange(cellid, run, status)

## create directories for each cell type
metadata_scharer %>%
    pull(cellid) %>%
    unique() %>%
    file.path("./results/scharer", .) %>%
    walk(~dir.create(., recursive = TRUE))

## junction files
juncfiles_scharer <- 
    "/lab-share/IM-Gutierrez-e2/Public/vitor/lupus/splicing/mapping/scharer/%s.junc" %>%
    sprintf(metadata_scharer$run)

if (all(file.exists(juncfiles_scharer))) {
    tibble(f = juncfiles_scharer, run = str_extract(juncfiles_scharer, "SRR\\d+")) %>%
	left_join(metadata_scharer) %>%
	group_by(cellid) %>%
	nest() %>%
	ungroup() %>%
	{walk2(.$cellid, .$data, 
	       ~write_lines(.y$f, file.path("./results/scharer", .x, "juncfiles.txt")))}
}

## groups file
metadata_scharer %>%
    group_by(cellid) %>%
    nest() %>%
    {walk2(.$cellid, .$data, 
           ~write_tsv(.y, file.path("./results/scharer", .x, "groups_file.txt"), col_names = FALSE))}

# Barnas
metadata_barnas <- 
    "/lab-share/IM-Gutierrez-e2/Public/External_datasets/Barnas/SraRunTable.txt" %>%
    read_csv() %>%
    mutate(cellid = case_when(flow_sort_selection == "CD19+ CD3- CD20+ IgD- CD27-" ~ "DN",
			      flow_sort_selection == "CD19+ CD3- CD20+ IgD+ CD27-" ~ "Naive",
			      flow_sort_selection == "CD19- CD14+ CD3-" ~ "Monocytes")) %>%
    filter(cellid != "Monocytes") %>%
    select(cellid, run = Run, status = disease_state) 

metadata_barnas %>%
    pull(cellid) %>%
    unique() %>%
    file.path("./results/barnas", .) %>%
    walk(~dir.create(., recursive = TRUE))

##
juncfiles_barnas <- 
    "/lab-share/IM-Gutierrez-e2/Public/vitor/lupus/splicing/mapping/barnas/%s.junc" %>%
    sprintf(metadata_barnas$run)

if ( all(file.exists(juncfiles_barnas)) ) {
    tibble(f = juncfiles_barnas, run = str_extract(juncfiles_barnas, "SRR\\d+")) %>%
	left_join(metadata_barnas) %>%
	group_by(cellid) %>%
	nest() %>%
	ungroup() %>%
	{walk2(.$cellid, .$data, 
	       ~write_lines(.y$f, file.path("./results/barnas", .x, "juncfiles.txt")))}
}


## groups file
metadata_barnas %>%
    arrange(cellid, status, run) %>%
    group_by(cellid) %>%
    nest() %>%
    {walk2(.$cellid, .$data, 
           ~write_tsv(.y, file.path("./results/barnas", .x, "groups_file.txt"), col_names = FALSE))}
   
# Andreoletti
samples_and <- 
    "/lab-share/IM-Gutierrez-e2/Public/External_datasets/Andreoletti/SRR_Acc_List_Bcells.txt" %>%
    read_lines()

juncfiles_and <- 
    "/lab-share/IM-Gutierrez-e2/Public/vitor/lupus/splicing/mapping/andreoletti/%s.junc" %>%
    sprintf(samples_and)

dir.create("./results/andreoletti/B", recursive = TRUE)

if ( all(file.exists(juncfiles_and)) ) {

    write_lines(juncfiles_and, "./results/andreoletti/B/juncfiles.txt") 
}

## groups file
read_tsv("../andreoletti/clusters_IFNg.tsv") %>%
    mutate(cluster = recode(cluster, "1" = "H", "2" = "L"),
	   cluster = factor(cluster, levels = c("L", "H"))) %>%
    arrange(cluster, sampleid) %>%
    write_tsv("./results/andreoletti/B/groups_file.txt", col_names = FALSE)

