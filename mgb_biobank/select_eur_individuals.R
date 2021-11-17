library(tidyverse)

# Metadata for 1000G samples
index_1000G <- "https://raw.githubusercontent.com/igsr/1000Genomes_data_indexes/master/data_collections/1000_genomes_project/1000genomes.sequence.index"

sample_annotation <- read_tsv(index_1000G, comment = "##") %>%
    select(sample_id = SAMPLE_NAME,
           population = POPULATION) %>%
    distinct()

# read PCA results into R
pca_genos <-
    "/lab-share/IM-Gutierrez-e2/Public/vitor/ase/mgb_biobank/results/plink_pca.eigenvec" %>%
    read_table(col_names = FALSE) %>%
    select(-1) %>%
    select(X2:X5) %>%
    setNames(c("sample_id", paste0("PC", 1:3)))

pca_eur <- pca_genos %>%
    inner_join(sample_annotation) %>%
    filter(population %in% c("CEU", "GBR", "FIN", "TSI", "IBS")) %>%
    select(sample_id, population, PC1:PC3)

eur_mu <- pca_eur %>%
    summarise_at(vars(PC1:PC3), mean)

pca_mgb <- pca_genos %>%
    anti_join(sample_annotation) %>%
    mutate(population = "MGB_biobank") %>%
    select(sample_id, population, PC1:PC3)

mgb_selected <- pca_mgb %>%
    mutate(z_pc1 = (PC1 - eur_mu$PC1)/sd(PC1),
       z_pc2 = (PC2 - eur_mu$PC2)/sd(PC2),
       z_pc3 = (PC3 - eur_mu$PC3)/sd(PC3)) %>%
    filter(between(z_pc1, -1, 1),
	   between(z_pc2, -1, 1),
	   between(z_pc3, -1, 1))

mgb_selected %>%
    pull(sample_id) %>%
    write_lines("./mgb_mosteuropean.txt")
