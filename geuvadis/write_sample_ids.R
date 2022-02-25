library(tidyverse)

samples_phase3 <- 
    "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel" %>%
    read_tsv(skip = 1, col_names = FALSE) %>%
    select(subject = X1, pop = X2, continent = X3, sex = X4)

geuvadis_info <- 
    "http://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/E-GEUV-1.sdrf.txt" %>%
    read_tsv() %>%
    select(name = `Source Name`, 
	   ena_id = `Comment[ENA_RUN]`) %>% 
    distinct() %>%
    inner_join(samples_phase3, by = c("name" = "subject"))

write_lines(geuvadis_info$ena_id, "./geuvadis_samples.txt")
