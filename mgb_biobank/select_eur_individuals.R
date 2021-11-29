library(tidyverse)

# Metadata for 1000G samples
index_1000G <- 
    file.path("https://raw.githubusercontent.com/igsr/1000Genomes_data_indexes",
	      "master/data_collections/1000_genomes_project",
	      "1000genomes.sequence.index")

sample_annotation <- read_tsv(index_1000G, comment = "##") %>%
    select(sample_id = SAMPLE_NAME,
           population = POPULATION) %>%
    distinct()

# read PCA results into R
pca_genos <-
    file.path("/lab-share/IM-Gutierrez-e2/Public/vitor",
	      "ase/mgb_biobank/results/allchr.merged.pruned.pca.eigenvec") %>%
    read_table(col_names = FALSE) %>%
    select(-1) %>%
    select(X2:X4) %>%
    setNames(c("sample_id", paste0("PC", 1:2)))

pca_eur <- pca_genos %>%
    inner_join(sample_annotation) %>%
    filter(population %in% c("CEU", "GBR", "FIN", "TSI", "IBS")) %>%
    select(sample_id, population, PC1:PC2)

eur_mu <- pca_eur %>%
    summarise_at(vars(PC1:PC2), list(mean = mean, sd = sd))

pca_mgb <- pca_genos %>%
    anti_join(sample_annotation) %>%
    mutate(population = "MGB_biobank") %>%
    select(sample_id, population, PC1:PC2)

mgb_selected <- pca_mgb %>%
    mutate(z_pc1 = abs(PC1 - eur_mu$PC1_mean)/eur_mu$PC1_sd,
	   z_pc2 = abs(PC2 - eur_mu$PC2_mean)/eur_mu$PC2_sd) %>%
    filter(between(z_pc1, 0, 3),
	   between(z_pc2, 0, 3))

mgb_selected %>%
    pull(sample_id) %>%
    write_lines("./results/mgb_eur.txt")
