library(tidyverse)

labshr <- "/lab-share/IM-Gutierrez-e2/Public/External_datasets" 
datasets <- c("Andreoletti", "Barnas", "Scharer")  

# Check sequence depth and remove samples with very low 
mapping_stats <- 
    file.path(labshr, 
	      c("Andreoletti/qc/multiqc_data/multiqc_fastqc.txt",
		"Barnas/qc/multiqc_data/multiqc_fastqc.txt",
		"Scharer/RNAseq/qc/multiqc_data/multiqc_fastqc.txt")) %>%
    setNames(datasets) %>%
    map_df(~read_tsv(.) %>%
               select(Sample, "Total Sequences", total_deduplicated_percentage) %>%
               mutate(unique = `Total Sequences` * total_deduplicated_percentage / 100,
                      duplicate = `Total Sequences` - unique) %>%
               select(Sample, unique, duplicate), .id = "dataset") %>%
    pivot_longer(unique:duplicate, names_to = "type", values_to = "n")

samples_rm1 <- mapping_stats %>%
    group_by(dataset, Sample) %>%
    summarise(n = sum(n)) %>%
    ungroup() %>%
    mutate(Sample = sub("^([^_]+).*$", "\\1", Sample)) %>%
    filter(n < 2e6) %>%
    distinct(dataset, sampleid = Sample)
    
# Remove samples that fail at many stats
samples_rm2 <- 
    file.path(labshr, 
	      c("Andreoletti/qc/multiqc_data/multiqc_fastqc.txt",
		"Barnas/qc/multiqc_data/multiqc_fastqc.txt",
		"Scharer/RNAseq/qc/multiqc_data/multiqc_fastqc.txt")) %>%
    map_df(~read_tsv(., col_types = c(.default = "c"))) %>% 
    select(Sample, basic_statistics:adapter_content) %>%
    pivot_longer(-Sample, names_to = "stat") %>%
    mutate(Sample = sub("^(SRR\\d+).*$", "\\1", Sample)) %>%
    group_by(Sample, stat) %>%
    summarise(value = ifelse(any(value == "fail"), "fail", "pass")) %>%
    ungroup() %>%
    count(Sample, value) %>%
    filter(value == "fail", n > 4)

samples_rm <- c(samples_rm1$sampleid, samples_rm2$Sample)

# Scharer
metadata_scharer <- 
    file.path(labshr, "Scharer/RNAseq/metadata/file_description.tsv") %>%
    read_tsv() %>%
    separate(id, c("status", "dummy", "cellid"), sep = "\\.") %>%
    select(cellid, run, status) %>%
    arrange(cellid, run, status) %>%
    filter(!run %in% samples_rm)

## create directories for each cell type
metadata_scharer %>%
    pull(cellid) %>%
    unique() %>%
    file.path("./results/scharer", .) %>%
    walk(~dir.create(., recursive = TRUE))

## junction files
juncfiles_scharer <- 
    "../../read_mapping/scharer/mapping/%s.junc" %>%
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
    file.path(labshr, "Barnas/SraRunTable.txt") %>%
    read_csv() %>%
    mutate(cellid = case_when(flow_sort_selection == "CD19+ CD3- CD20+ IgD- CD27-" ~ "DN",
			      flow_sort_selection == "CD19+ CD3- CD20+ IgD+ CD27-" ~ "Naive",
			      flow_sort_selection == "CD19- CD14+ CD3-" ~ "Monocytes")) %>%
    filter(cellid != "Monocytes") %>%
    select(cellid, run = Run, status = disease_state) %>%
    filter(!run %in% samples_rm)

metadata_barnas %>%
    pull(cellid) %>%
    unique() %>%
    file.path("./results/barnas", .) %>%
    walk(~dir.create(., recursive = TRUE))

##
juncfiles_barnas <- 
    "../../read_mapping/barnas/mapping/%s.junc" %>%
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
samples_andre <- 
    file.path(labshr, "Andreoletti/SRR_Acc_List_Bcells.txt") %>%
    read_lines() %>%
    tibble(run = .) %>%
    filter(!run %in% samples_rm)

juncfiles_andre <- 
    "../../read_mapping/andreoletti/mapping/%s.junc" %>%
    sprintf(samples_andre$run)

dir.create("./results/andreoletti/B", recursive = TRUE)

if ( all(file.exists(juncfiles_andre)) ) {

    write_lines(juncfiles_andre, "./results/andreoletti/B/juncfiles.txt") 
}

## groups file
ifn_andre <- read_tsv("../../read_mapping/andreoletti/clusters_IFNg.tsv") %>%
    mutate(cluster = recode(cluster, "1" = "H", "2" = "L"),
	   cluster = factor(cluster, levels = c("L", "H"))) %>%
    arrange(cluster, sampleid) %>%
    filter(!sampleid %in% samples_rm) 

write_tsv(ifn_andre, "./results/andreoletti/B/groups_file.txt", col_names = FALSE)

metadata_andre <- left_join(samples_andre, ifn_andre, by = c("run" = "sampleid")) %>%
    mutate(cellid = "B") %>%
    select(cellid, run, status = cluster)

# Create file for slurm batch jobs
bind_rows(
    scharer = metadata_scharer, 
    barnas = metadata_barnas,
    andreoletti = metadata_andre,
    .id = "dataset") %>%
    distinct(dataset, cellid) %>%
    write_tsv("./slurm_info.txt", col_names = FALSE)
