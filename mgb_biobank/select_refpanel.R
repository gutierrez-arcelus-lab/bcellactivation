library(tidyverse)

# Metadata for 1000G samples
index_1000G <- 
    "https://raw.githubusercontent.com/igsr/1000Genomes_data_indexes/master/data_collections/1000_genomes_project/1000genomes.sequence.index"

sample_annotation <- read_tsv(index_1000G, comment = "##") %>%
    select(sample_id = SAMPLE_NAME,
           population = POPULATION) %>%
    distinct()

sample_annotation %>%
    #filter(population %in% c("CEU", "YRI", "CHS", "ITU")) %>%
    filter(population %in% c("CEU", "YRI", "CHS")) %>%
    arrange(sample_id) %>%
    mutate(id = paste0(sample_id, "_", sample_id)) %>%
    pull(id) %>%
    write_lines("./refpanel_ids.txt")
