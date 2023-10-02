library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(uwot)

batches <- sprintf("04%02d", c(1:8, 10))

## PCA on genotype data

### Metadata for 1000G samples
index_1000G <- 
    "https://raw.githubusercontent.com/igsr/1000Genomes_data_indexes/master/data_collections/1000_genomes_project/1000genomes.sequence.index"

sample_annotation <- read_tsv(index_1000G, comment = "##") %>%
    select(sample_id = SAMPLE_NAME,
           population = POPULATION) %>%
    distinct()


### Colors
afr_cols <- brewer.pal("Oranges", n = 9)[3:9] %>%
    setNames(c("ACB", "ESN", "GWD", "LWK", "MSL", "YRI", "ASW"))

eur_cols <- brewer.pal("Blues", n = 9)[5:9] %>%
    setNames(c("GBR", "IBS", "TSI", "CEU", "FIN"))

sas_cols <- brewer.pal("Greens", n = 9)[5:9] %>%
    setNames(c("BEB", "PJL", "GIH", "ITU", "STU"))

eas_cols <- brewer.pal("Purples", n = 9)[5:9] %>%
    setNames(c("CHB", "CHS", "CDX", "KHV", "JPT"))

amr_cols <- c("lightpink1", "hotpink", "hotpink3", "deeppink") %>%
    setNames(c("MXL", "CLM", "PEL", "PUR"))

mgb_cols <- c("grey", "turquoise", "pink", "brown") %>%
    setNames(c("MGB_biobank", "MGB_eur", "MGB_amr", "MGB_aa"))

all_cols <- c(mgb_cols, afr_cols, eur_cols, sas_cols, eas_cols, amr_cols)

##### Pca for genotyped SNPs only
pca_genos <- 
    "./results/VCF/genotyped/allchr.isec.merged.pruned.pca.eigenvec" %>%
    read_table(col_names = FALSE) %>%
    setNames(c("id1", "id2", paste0("PC", 1:22)))

pca_kgp <- pca_genos %>%
    inner_join(sample_annotation, by = c("id2" = "sample_id")) %>%
    select(sample_id = id2, population, starts_with("PC"))

pca_mgb <- pca_genos %>%
    anti_join(sample_annotation, by = c("id2" = "sample_id")) %>%
    mutate(population = "MGB_biobank") %>%
    unite("sample_id", c("id1", "id2"), sep = "_") %>%
    select(sample_id, population, starts_with("PC"))


pca_df <- bind_rows(pca_mgb, pca_kgp) %>%
    mutate(population = factor(population, levels = names(all_cols))) %>%
    mutate(dataset = ifelse(grepl("MGB", population), "MGB", "1000 Genomes"),
           dataset = factor(dataset, levels = rev(c("1000 Genomes", "MGB")))) 

## Kmeans
add_cluster <- function(k, pca_data) {
    
    clusters <- kmeans(pca_matrix, k)$cluster
    
    pca_data %>%
        add_column(cluster = clusters)
}

pca_matrix <- pca_df %>%
    select(PC1:PC6) %>%
    as.matrix()


clusters_df <- map_df(1:10, ~add_cluster(k = ., pca_data = pca_df), .id = "k") %>%
    mutate(k = as.integer(k)) %>%
    arrange(k)

clusters_continent <- clusters_df %>%
    filter(dataset != "MGB") %>%
    mutate(continent = case_when(population %in% names(afr_cols) ~ "AFR",
                                 population %in% names(eur_cols) ~ "EUR",
                                 population %in% names(sas_cols) ~ "SAS",
                                 population %in% names(eas_cols) ~ "EAS",
                                 population %in% names(amr_cols) ~ "AMR"),
           continent = factor(continent,
                              levels = c("AFR", "EUR", "SAS", "EAS", "AMR"))) %>%
    group_by(k, cluster) %>%
    count(continent) %>%
    mutate(p = n/sum(n)) %>%
    slice(which.max(n)) %>%
    select(k, cluster, continent, p) %>%
    ungroup()

clusters_continent %>%
    group_by(k) %>%
    summarise(p = mean(p)) %>%
    ggplot(aes(factor(k), p)) +
    geom_point() +
    geom_line(aes(group = 1)) +
    theme_bw() +
    labs(x = "K", y = "Precision")

continent_colors <- all_cols[c("YRI", "TSI", "PJL", "CDX", "PEL")] %>%
    setNames(levels(clusters_continent$continent)) %>%
    c("NA" = "grey")

clusters_plot_df <- clusters_df %>%
    left_join(select(clusters_continent, -p), by = c("k", "cluster")) %>%
    filter(dataset == "MGB") %>%
    mutate(continent = as.character(continent),
           continent = factor(continent, levels = names(continent_colors)))

kmeans_df <- pca_df %>% 
    filter(dataset == "1000 Genomes") %>%
    mutate(continent = case_when(population %in% names(afr_cols) ~ "AFR",
                                 population %in% names(eur_cols) ~ "EUR",
                                 population %in% names(sas_cols) ~ "SAS",
                                 population %in% names(eas_cols) ~ "EAS",
                                 population %in% names(amr_cols) ~ "AMR"),
           continent = factor(continent,
                              levels = c("AFR", "EUR", "SAS", "EAS", "AMR"))) %>% 
    bind_rows(clusters_plot_df %>% filter(dataset == "MGB", k == 10) %>% select(-k, -cluster))

pca_kmeans_plot_1 <- kmeans_df %>%
    ggplot(aes(PC1, PC2, color = continent, alpha = dataset)) +
    geom_point(size = .5) +
    scale_color_manual(values = continent_colors) +
    scale_alpha_manual(values = c("MGB" = .5, "1000 Genomes" = 1)) +
    facet_wrap(~dataset) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.text = element_text(face = "bold", size = 11),
          panel.spacing = unit(2, "lines")) +
    guides(color = guide_legend(override.aes = list(size = 2), ncol = 1),
           alpha = "none")

pca_kmeans_plot_2 <- kmeans_df %>%
    ggplot(aes(PC3, PC4, color = continent, alpha = dataset)) +
    geom_point(size = .5) +
    scale_color_manual(values = continent_colors) +
    scale_alpha_manual(values = c("MGB" = .5, "1000 Genomes" = 1)) +
    facet_wrap(~dataset) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.text = element_text(face = "bold", size = 11),
          panel.spacing = unit(2, "lines")) +
    guides(color = guide_legend(override.aes = list(size = 2), ncol = 1),
           alpha = "none")

pca_kmeans_plot_3 <- kmeans_df %>%
    ggplot(aes(PC5, PC6, color = continent, alpha = dataset)) +
    geom_point(size = .5) +
    scale_color_manual(values = continent_colors) +
    scale_alpha_manual(values = c("MGB" = .5, "1000 Genomes" = 1)) +
    facet_wrap(~dataset) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.text = element_text(face = "bold", size = 11),
          panel.spacing = unit(2, "lines")) +
    guides(color = guide_legend(override.aes = list(size = 2), ncol = 1),
           alpha = "none")


pca_kmeans_grid <- plot_grid(pca_kmeans_plot_1 + theme(legend.position = "none"),
                             pca_kmeans_plot_2 + theme(legend.position = "none"),
                             pca_kmeans_plot_3 + theme(legend.position = "none"),
                             ncol = 1)

plot_grid(pca_kmeans_grid, get_legend(pca_kmeans_plot_1),
          nrow = 1, rel_widths = c(1, .3))

ggsave("./plots/pca_kmeans.png", width = 8, height = 7.55)

kmeans_df %>%
    filter(dataset == "MGB") %>%
    count(continent)

## UMAP
"./results/VCF/genotyped/allchr.isec.merged.pruned.pca.eigenval" %>%
    read_lines() %>%
    tibble(eigenval = as.numeric(.)) %>%
    select(eigenval) %>%
    mutate(var_exp = eigenval/sum(eigenval) * 100)

mgb_umap <- select(pca_df, PC1:PC6) %>%
    umap()

umap_df <- bind_cols(select(pca_df, sample_id, dataset, population), as.data.frame(mgb_umap))

plot_umap_pc6 <- umap_df %>%
    ggplot(aes(V1, V2, color = population)) +
    geom_point(size = .5) +
    scale_color_manual(values = all_cols[-(2:4)]) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 2), ncol = 2),
           alpha = "none") +
    labs(x = "UMAP 1", y = "UMAP 2")

ggsave("./plots/umap_6pcs.png", plot_umap_pc6, width = 8, height = 6)


