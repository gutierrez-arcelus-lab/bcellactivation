library(tidyverse)
library(readxl)

samples_phase3 <- 
    "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel" %>%
    read_tsv(skip = 1, col_names = FALSE) %>%
    select(subject = X1, pop = X2, continent = X3, sex = X4)

geuvadis_info <- 
    "http://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/E-GEUV-1.sdrf.txt" %>%
    read_tsv() %>% 
    select(name = `Source Name`, 
	   ena_id = `Comment[ENA_RUN]`,
	   lab = Performer) %>% 
    distinct() %>%
    inner_join(samples_phase3, by = c("name" = "subject"))

write_lines(geuvadis_info$ena_id, "./geuvadis_samples.txt")

geuvadis_info %>%
    filter(continent == "EUR") %>%
    pull(name) %>%
    write_lines("./geuvadis_1000G_ids.txt")


ebv <- read_excel("./ebv_copynumbers_pone.0179446.xlsx") %>%
    select(-pop)

geuvadis_info %>%
    filter(continent == "EUR") %>%
    left_join(ebv, by = c("name" = "samples")) %>%
    select(sampleid = ena_id, pop, sex, lab, ebv_load = `EBV load`) %>%
    write_tsv("./geuvadis_covariates.tsv")

geuvadis_info %>%
    select(ena_id, name) %>%
    write_tsv("./geuvadis_ids_conversionkey.tsv")
