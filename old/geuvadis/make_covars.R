library(tidyverse)

pca <- read_delim("./VCF/allchr.1000G.pca", delim = " ") %>%
    rename(id = SampleID) %>%
    mutate(id = str_extract(id, "PC\\d+$")) %>%
    filter(id %in% c("PC1", "PC2"))

covars_df <- read_tsv("./geuvadis_covariates.tsv")

lab_df <- covars_df %>%
    select(kgp_id, lab) %>%
    mutate(lab = paste0("lab_", lab),
           i = 1) %>%
    complete(kgp_id, lab, fill = list(i = 0)) %>%
    pivot_wider(names_from = lab, values_from = i)

sex_df <- covars_df %>%
    mutate(sex = ifelse(sex == "female", 0, 1)) %>%
    select(kgp_id, sex)

ebv_df <- covars_df %>%
    select(kgp_id, ebv_load)

covar_matrix <- left_join(lab_df, sex_df) %>%
    left_join(ebv_df) %>%
    arrange(kgp_id) %>%
    pivot_longer(-kgp_id, names_to = "id") %>%
    pivot_wider(names_from = kgp_id, values_from = value) %>%
    bind_rows(pca[names(.)], .)

write_tsv(covar_matrix, "./qtltools_cov.txt")
