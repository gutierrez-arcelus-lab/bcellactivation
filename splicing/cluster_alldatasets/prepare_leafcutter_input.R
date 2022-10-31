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
	      c(
		"Andreoletti/qc/multiqc_data/multiqc_fastqc.txt",
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

# Metadata for the datasets
## Scharer
metadata_scharer <- 
    file.path(labshr, "Scharer/RNAseq/metadata/file_description.tsv") %>%
    read_tsv() %>%
    separate(id, c("status", "dummy", "cellid"), sep = "\\.") %>%
    select(cellid, run, status) %>%
    arrange(cellid, run, status) %>%
    filter(!run %in% samples_rm)

## Barnas
metadata_barnas <- 
    file.path(labshr, "Barnas/SraRunTable.txt") %>%
    read_csv() %>%
    mutate(cellid = case_when(flow_sort_selection == "CD19+ CD3- CD20+ IgD- CD27-" ~ "DN",
			      flow_sort_selection == "CD19+ CD3- CD20+ IgD+ CD27-" ~ "Naive",
			      flow_sort_selection == "CD19- CD14+ CD3-" ~ "Monocytes")) %>%
    filter(cellid != "Monocytes") %>%
    select(cellid, run = Run, status = disease_state) %>%
    filter(!run %in% samples_rm) 

## Andreoletti
metadata_andre <- "../../read_mapping/andreoletti/clusters_IFNg.tsv" %>%
    read_tsv() %>%
    mutate(cluster = recode(cluster, "1" = "H", "2" = "L"),
           cluster = factor(cluster, levels = c("L", "H"))) %>%
    arrange(cluster, sampleid) %>%
    add_column(cellid = "B", .before = 1) %>%
    select(cellid, run = sampleid, status = cluster) %>%
    filter(!run %in% samples_rm)


# Create results directories for each dataset and cell type
metadata_scharer %>%
    pull(cellid) %>%
    unique() %>%
    file.path("./results/scharer", .) %>%
    walk(~dir.create(., recursive = TRUE))

metadata_barnas %>%
    pull(cellid) %>%
    unique() %>%
    file.path("./results/barnas", .) %>%
    walk(~dir.create(., recursive = TRUE))

metadata_andre %>%
    pull(cellid) %>%
    unique() %>%
    file.path("./results/andreoletti", .) %>%
    walk(~dir.create(., recursive = TRUE))


# Write groups file defining which samples belong to each group
metadata_scharer %>%
    group_by(cellid) %>%
    nest() %>%
    {walk2(.$cellid, .$data, 
           ~write_tsv(.y, file.path("./results/scharer", .x, "groups_file.txt"), col_names = FALSE))}

metadata_barnas %>%
    arrange(cellid, status, run) %>%
    group_by(cellid) %>%
    nest() %>%
    {walk2(.$cellid, .$data, 
           ~write_tsv(.y, file.path("./results/barnas", .x, "groups_file.txt"), col_names = FALSE))}

metadata_andre %>%
    arrange(cellid, status, run) %>%
    group_by(cellid) %>%
    nest() %>%
    {walk2(.$cellid, .$data, 
           ~write_tsv(.y, file.path("./results/andreoletti", .x, "groups_file.txt"), col_names = FALSE))}

# Create file for slurm batch jobs
bind_rows(
    scharer = metadata_scharer, 
    barnas = metadata_barnas,
    andreoletti = metadata_andre,
    .id = "dataset") %>%
    distinct(dataset, cellid) %>%
    write_tsv("./slurm_info.txt", col_names = FALSE) 


# Junction files
juncfiles_scharer <- 
    "../read_mapping/scharer/mapping/%s.junc" %>%
    sprintf(metadata_scharer$run)

juncfiles_barnas <- 
    "../read_mapping/barnas/mapping/%s.junc" %>%
    sprintf(metadata_barnas$run)

juncfiles_andre <- 
    "../read_mapping/andreoletti/mapping/%s.junc" %>%
    sprintf(metadata_andre$run)

junctionfiles <- c(juncfiles_scharer, juncfiles_barnas, juncfiles_andre)

if (all(file.exists(junctionfiles)) ) {
    write_lines(junctionfiles, "./juncfiles.txt")
    cat("Done!\n")
} else {
    stop("missing junction files")
}
