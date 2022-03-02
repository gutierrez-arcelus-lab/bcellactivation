library(tidyverse)
library(peer)


covars_df <- read_tsv("./geuvadis_covariates.tsv") %>%
    filter_all(~!is.na(.))

pop_df <- covars_df %>% 
    select(sampleid, pop) %>%
    mutate(pop = paste0("pop_", pop),
	   i = 1) %>%
    complete(sampleid, pop, fill = list(i = 0)) %>%
    pivot_wider(names_from = pop, values_from = i)

lab_df <- covars_df %>%
    select(sampleid, lab) %>%
    mutate(lab = paste0("lab_", lab),
	   i = 1) %>%
    complete(sampleid, lab, fill = list(i = 0)) %>%
    pivot_wider(names_from = lab, values_from = i)

sex_df <- covars_df %>%
    mutate(sex = ifelse(sex == "female", 0, 1)) %>%
    select(sampleid, sex)

ebv_df <- covars_df %>%
    select(sampleid, ebv_load)

covar_matrix <- left_join(pop_df, lab_df) %>%
    left_join(sex_df) %>%
    left_join(ebv_df) %>%
    column_to_rownames("sampleid") %>%
    data.matrix()

phenotypes <- 
    read_tsv("./geuvadis_salmon_quants_ebv.bed") %>%
    select(id, starts_with("ERR")) %>%
    pivot_longer(-id, names_to = "sampleid") %>%
    filter(sampleid %in% covars_df$sampleid) %>%
    pivot_wider(names_from = id) %>%
    column_to_rownames("sampleid") %>%
    data.matrix()

model <- PEER()
PEER_setPhenoMean(model, phenotypes)
PEER_setCovariates(model, covar_matrix)
PEER_setNk(model, 50)
PEER_setNmax_iterations(model, 1e6)
PEER_setAdd_mean(model, TRUE)
PEER_update(model)

write_rds(model, "./results/peer_model.rds")
