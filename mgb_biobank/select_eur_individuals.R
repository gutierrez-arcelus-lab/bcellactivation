library(tidyverse)

# function to extract MGB individuals falling in up to 1 sd from 
# the average of any European population
filter_pca <- function(dat, pc) {
    pcmin <- paste0(pc, "_min") 
    pcmax <- paste0(pc, "_max")

    dat %>%
    filter(between(PC1, eur_summary[[pcmin]][1], eur_summary[[pcmax]][1]) |
	   between(PC1, eur_summary[[pcmin]][2], eur_summary[[pcmax]][2]) |
	   between(PC1, eur_summary[[pcmin]][3], eur_summary[[pcmax]][3]) |
	   between(PC1, eur_summary[[pcmin]][4], eur_summary[[pcmax]][4]) |
	   between(PC1, eur_summary[[pcmin]][5], eur_summary[[pcmax]][5]))
}

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
	      "ase/mgb_biobank/results/plink_pca.eigenvec") %>%
    read_table(col_names = FALSE) %>%
    select(-1) %>%
    select(X2:X5) %>%
    setNames(c("sample_id", paste0("PC", 1:3)))

pca_eur <- pca_genos %>%
    inner_join(sample_annotation) %>%
    filter(population %in% c("CEU", "GBR", "FIN", "TSI", "IBS")) %>%
    select(sample_id, population, PC1:PC3)

#eur_summary <- pca_eur %>%
#    group_by(population) %>%
#    summarise_at(vars(PC1:PC3), 
#		 list(min = function(x) mean(x) - (2 * sd(x)), 
#		      max = function(x) mean(x) + (2 * sd(x))))
#

eur_mu <- pca_eur %>%
    summarise_at(vars(PC1:PC3), mean)

pca_mgb <- pca_genos %>%
    anti_join(sample_annotation) %>%
    mutate(population = "MGB_biobank") %>%
    select(sample_id, population, PC1:PC3)

#mgb_selected <- pca_mgb %>%
#    filter_pca("PC1") %>%
#    filter_pca("PC2") %>%
#    filter_pca("PC3")
#
mgb_selected <- pca_mgb %>%
    mutate(z_pc1 = abs(PC1 - eur_mu$PC1)/sd(PC1),
	   z_pc2 = abs(PC2 - eur_mu$PC2)/sd(PC2),
	   z_pc3 = abs(PC3 - eur_mu$PC3)/sd(PC3)) %>%
    filter(between(z_pc1, 0, 0.5),
	   between(z_pc2, 0, 0.5),
	   between(z_pc3, 0, 0.5))

mgb_selected %>%
    pull(sample_id) %>%
    write_lines("./mgb_mosteuropean.txt")

#bind_rows("mgb" = pca_mgb, "1000G" = pca_eur, .id = "dataset") %>%
#    mutate(most_eur = ifelse(sample_id %in% mgb_selected$sample_id, 1L, 0L)) %>%
#    write_tsv("./results/mgb_pcavalues.tsv")
