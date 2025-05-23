library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(ggrepel)
library(ggthemes)
library(ggsci)
library(ggbeeswarm)

## Select females
batches <- sprintf("04%02d", c(1:8, 10))

het <- paste0("/temp_work/ch229163/VCF/chrX.", batches, ".genotyped.het") %>%
    setNames(batches) %>%
    map_dfr(read_tsv, .id = "batch") %>%
    mutate(hom = `O(HOM)`/N_SITES) 

sex_1 <- ggplot(het, aes(hom)) +
    geom_histogram(bins = 50) +
    geom_vline(xintercept = .95, linetype = 2) +
    scale_x_continuous(labels = function(x) round(x, 2)) +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank()) +
    labs(x = "Homozygosity")

sex_2 <- het %>%
    mutate(sex = ifelse(hom > .95, "male", "female")) %>%
    count(sex) %>%
    mutate(p = n/sum(n)) %>%
    ggplot(aes(x = 1, y = p, fill = sex)) +
    geom_col(position = "fill") +
    geom_hline(yintercept = .5, linetype = 2) +
    scale_fill_npg(labels = c("< 0.985", ">0.985")) +
    scale_y_continuous(breaks = c(0, .5, 1), labels = scales::percent) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_blank(),
          legend.title = element_text(size = 9)) +
    labs(x = NULL, y = NULL, fill = "Homozygosity") +
    guides(color = guide_legend(override.aes = list(size = 2)))

plot_grid(sex_1, NULL, sex_2, nrow = 1, rel_widths = c(1, .1, .6))

ggsave("./plots/chrX_het.png", width = 5)

# Write ID files
females_df <- het %>%
    filter(hom < .95) %>%
    select(1:2) 

females_df %>%
    group_nest(batch) %>%
    mutate(out = paste0("./results/females_", batch, ".txt"),
           data = map(data, unlist)) %>%
    walk2(.x = .$data, .y = .$out, .f = ~write_lines(.x, .y))

##

females_df <- paste0("./results/females_", batches, ".txt") %>%
    setNames(batches) %>%
    map_df(~tibble(sample_id = read_lines(.)), .id = "batch")

##

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


### read PCA results into R
pca_genos <- "./results/VCF/allchr.isec.merged.pruned.pca.eigenvec" %>%
    read_table(col_names = FALSE) %>%
    select(X1:X6) %>%
    setNames(c("id1", "id2", paste0("PC", 1:4)))

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

pca_thresholds <- pca_kgp %>%
    filter(population %in% c("GBR", "CEU", "TSI", "IBS", "FIN")) %>%
    pivot_longer(PC1:PC4, names_to = "pc") %>%
    group_by(pc) %>%
    slice(-which.max(value)) %>%
    slice(-which.min(value)) %>%
    summarise(value = range(value)) %>%
    mutate(i = c("min", "max"),
           value = value * 1.1) %>%
    ungroup() %>%
    pivot_wider(names_from = pc, values_from = value) %>%
    mutate(dataset = "MGB",
           dataset = factor(dataset, levels = c("1000 Genomes", "MGB")))


pca_plot_1 <- pca_df %>%
    ggplot(aes(PC1, PC2, color = population, alpha = dataset)) +
    geom_point(size = .5) +
    scale_color_manual(values = all_cols) +
    scale_alpha_manual(values = c("MGB" = .1, "1000 Genomes" = 1)) +
    facet_wrap(~dataset) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.text = element_text(face = "bold", size = 11),
          panel.spacing = unit(2, "lines")) +
    guides(color = guide_legend(override.aes = list(size = 2), ncol = 2),
           alpha = FALSE)

pca_plot_2 <- pca_df %>%
    ggplot(aes(PC3, PC4, color = population, alpha = dataset)) +
    geom_point(size = .5) +
    scale_color_manual(values = all_cols) +
    scale_alpha_manual(values = c("MGB" = .1, "1000 Genomes" = 1)) +
    facet_wrap(~dataset) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.text = element_text(face = "bold", size = 11),
          panel.spacing = unit(2, "lines")) +
    guides(color = guide_legend(override.aes = list(size = 2), ncol = 2),
           alpha = FALSE)

pca_grid <- plot_grid(pca_plot_1 + theme(legend.position = "none"),
                      NULL,
                      pca_plot_2 + theme(legend.position = "none"),
                      ncol = 1, rel_heights = c(1, .05, 1))

pca_plot <- plot_grid(pca_grid, get_legend(pca_plot_1),
                      nrow = 1, rel_widths = c(1, .3))

ggsave("./plots/pca.png", pca_plot, width = 8, height = 5)


pca_df %>%
    select(sample_id, PC1, PC2) %>%
    pivot_longer(PC1:PC2, names_to = "PC") %>%
    ggplot(aes(value)) +
    geom_histogram(bins = 50) +
    facet_wrap(~PC, scales = "free") +
    theme_bw() +
    theme(panel.grid = element_blank())

ggsave("./plots/pca_histogram.png", width = 5, height = 2.5)


## Select based on PC1:PC4 coordinates
mgb_eur_inds <- pca_df %>%
    filter(dataset == "MGB") %>%
    filter(between(PC1, pca_thresholds$PC1[1], pca_thresholds$PC1[2]),
           between(PC2, pca_thresholds$PC2[1], pca_thresholds$PC2[2]),
           between(PC3, pca_thresholds$PC3[1], pca_thresholds$PC3[2]),
           between(PC4, pca_thresholds$PC4[1], pca_thresholds$PC4[2])) %>%
    select(sample_id)

inner_join(mgb_eur_inds, females_df, by = c("sample_id"))  %>%
    select(batch, sample_id) %>%
    group_nest(batch) %>%
    mutate(out = paste0("./results/eur_females_", batch, ".txt"),
           data = map(data, unlist)) %>%
    walk2(.x = .$data, .y = .$out, .f = ~write_lines(.x, .y))


pca_eur_df <- pca_df %>%
    select(sample_id, dataset, population, PC1:PC4) %>%
    mutate(population = as.character(population),
           population = ifelse(sample_id %in% mgb_eur_inds$sample_id, "MGB_eur", population),
           population = factor(population, levels = names(all_cols)))

all_cols[1:2] <- c("grey", "black")

pca_eur_plot_1 <- pca_eur_df %>%
    ggplot(aes(PC1, PC2, color = population, alpha = dataset)) +
    geom_point(size = .5) +
    geom_vline(data = pca_thresholds, 
               aes(xintercept = PC1), linetype = 2, alpha = .5, color = "red") +
    geom_hline(data = pca_thresholds, 
               aes(yintercept = PC2), linetype = 2, alpha = .5, color = "red") +
    scale_color_manual(values = all_cols) +
    scale_alpha_manual(values = c("MGB" = .1, "1000 Genomes" = 1)) +
    facet_wrap(~dataset) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.text = element_text(face = "bold", size = 11),
          panel.spacing = unit(2, "lines")) +
    guides(color = guide_legend(override.aes = list(size = 2), ncol = 2),
           alpha = FALSE)

pca_eur_plot_2 <- pca_eur_df %>%
    ggplot(aes(PC3, PC4, color = population, alpha = dataset)) +
    geom_point(size = .5) +
    geom_vline(data = pca_thresholds, 
               aes(xintercept = PC3), linetype = 2, alpha = .5, color = "red") +
    geom_hline(data = pca_thresholds, 
               aes(yintercept = PC4), linetype = 2, alpha = .5, color = "red") +
    scale_color_manual(values = all_cols) +
    scale_alpha_manual(values = c("MGB" = .1, "1000 Genomes" = 1)) +
    facet_wrap(~dataset) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.text = element_text(face = "bold", size = 11),
          panel.spacing = unit(2, "lines")) +
    guides(color = guide_legend(override.aes = list(size = 2), ncol = 2),
           alpha = FALSE)


pca_eur_grid <- plot_grid(pca_eur_plot_1 + theme(legend.position = "none"),
                          NULL,
                          pca_eur_plot_2 + theme(legend.position = "none"),
                          ncol = 1, rel_heights = c(1, .05, 1))

pca_eur_plot <- plot_grid(pca_eur_grid, get_legend(pca_eur_plot_1),
                          nrow = 1, rel_widths = c(1, .3))

ggsave("./plots/pca_eur.png", pca_eur_plot, width = 8, height = 5)


### Select admixed individuals
get_pop_pca <- function(pca_data, pop = c("eur", "aa", "amr")) {
    
    pop <- match.arg(pop)
    
    if (pop == "eur") pops <- c("CEU", "GBR", "TSI", "FIN")
    if (pop == "aa") pops <- c("ACB", "ESN", "GWD", "LWK", "MSL", "YRI", "ASW")
    if (pop == "amr") pops <- c("MXL", "CLM", "PEL", "PUR")

     pca_thresholds <- pca_data %>%
         filter(population %in% all_of(pops)) %>%
         pivot_longer(starts_with("PC"), names_to = "pc") %>%
         group_by(pc) %>%
         mutate(q = percent_rank(value)) %>%
         filter(q >= .01, q <= .99) %>%
         summarise(tibble(min = min(value), max = max(value))) %>%
         ungroup()
     
     mgb_inds <- pca_data %>%
         filter(dataset == "MGB") %>%
         pivot_longer(starts_with("PC"), names_to = "pc") %>%
         inner_join(pca_thresholds) %>%
         group_by(sample_id) %>%
         filter(all(value >= min & value <= max)) %>%
         ungroup() %>%
         distinct(sample_id)
    
    pca_data %>%
        inner_join(mgb_inds) %>% 
        mutate(population = paste0("MGB_", pop))
}

plot_pca <- function(pca_data) {
    
    pca_plot_1 <- pca_data %>%
        ggplot(aes(PC1, PC2, color = population, alpha = dataset)) +
        geom_point(size = .5) +
        scale_color_manual(values = all_cols) +
        scale_alpha_manual(values = c("MGB" = .5, "1000 Genomes" = 1)) +
        facet_wrap(~dataset) +
        theme_bw() +
        theme(panel.grid = element_blank(),
              strip.text = element_text(face = "bold", size = 11),
              panel.spacing = unit(2, "lines")) +
        guides(color = guide_legend(override.aes = list(size = 2), ncol = 2),
               alpha = "none")
    
    pca_plot_2 <- pca_data %>%
        ggplot(aes(PC3, PC4, color = population, alpha = dataset)) +
        geom_point(size = .5) +
        scale_color_manual(values = all_cols) +
        scale_alpha_manual(values = c("MGB" = .5, "1000 Genomes" = 1)) +
        facet_wrap(~dataset) +
        theme_bw() +
        theme(panel.grid = element_blank(),
              strip.text = element_text(face = "bold", size = 11),
              panel.spacing = unit(2, "lines")) +
        guides(color = guide_legend(override.aes = list(size = 2), ncol = 2),
               alpha = "none")
    
    
    pca_grid <- plot_grid(pca_plot_1 + theme(legend.position = "none"),
                          NULL,
                          pca_plot_2 + theme(legend.position = "none"),
                          ncol = 1, rel_heights = c(1, .05, 1))
    
    plot_grid(pca_grid, get_legend(pca_plot_1),
              nrow = 1, rel_widths = c(1, .3))
}

pca_aa <- get_pop_pca(pca_df, pop = "aa")
pca_amr <- get_pop_pca(pca_df, pop = "amr")

pca_df %>%
    anti_join(pca_aa, by = "sample_id") %>%
    anti_join(pca_amr, by = "sample_id") %>%
    bind_rows(pca_aa) %>%
    bind_rows(pca_amr) %>%
    plot_pca()

## Kmeans

pca_matrix <- typed_pca_df %>%
    select(starts_with("PC")) %>%
    as.matrix()

add_cluster <- function(k, pca_data) {
    
    clusters <- kmeans(pca_matrix, k)$cluster
    
    pca_data %>%
        add_column(cluster = clusters)
}

clusters_df <- map_df(1:10, ~add_cluster(k = ., pca_data = typed_pca_df), .id = "k") %>%
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
    labs(x = "K", y = "Homogeneity")

continent_colors <- all_cols[c("YRI", "TSI", "PJL", "CDX", "PEL")] %>%
    setNames(levels(clusters_continent$continent)) %>%
    c("NA" = "grey")

clusters_plot_df <- clusters_df %>%
    left_join(select(clusters_continent, -p), by = c("k", "cluster")) %>%
    filter(dataset == "MGB") %>%
    mutate(continent = as.character(continent),
           continent = factor(continent, levels = names(continent_colors)))

kmeans_df <- typed_pca_df %>% 
    filter(dataset == "1000 Genomes") %>%
    mutate(continent = case_when(population %in% names(afr_cols) ~ "AFR",
                                 population %in% names(eur_cols) ~ "EUR",
                                 population %in% names(sas_cols) ~ "SAS",
                                 population %in% names(eas_cols) ~ "EAS",
                                 population %in% names(amr_cols) ~ "AMR"),
           continent = factor(continent,
                              levels = c("AFR", "EUR", "SAS", "EAS", "AMR"))) %>% 
    bind_rows(clusters_plot_df %>% filter(dataset == "MGB", k == 9) %>% select(-k, -cluster))

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


pca_kmeans_grid <- plot_grid(pca_kmeans_plot_1 + theme(legend.position = "none"),
                             NULL,
                             pca_kmeans_plot_2 + theme(legend.position = "none"),
                             ncol = 1, rel_heights = c(1, .05, 1))

plot_grid(pca_kmeans_grid, get_legend(pca_kmeans_plot_1),
          nrow = 1, rel_widths = c(1, .3))

ggsave("./plots/pca_kmeans.png", width = 8, height = 5)

kmeans_df %>%
    filter(dataset == "MGB") %>%
    count(continent)

kmeans_df %>%
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

ggsave("./plots/pca_kmeans_pc5and6.png", width = 8, height = 2.5)


##### Pca for genotyped SNPs only
typed_pca_genos <- 
    "./results/VCF/genotyped/allchr.isec.merged.pruned.pca.eigenvec" %>%
    read_table(col_names = FALSE) %>%
    select(X1:X8) %>%
    setNames(c("id1", "id2", paste0("PC", 1:6)))

typed_pca_kgp <- typed_pca_genos %>%
    inner_join(sample_annotation, by = c("id2" = "sample_id")) %>%
    select(sample_id = id2, population, starts_with("PC"))

typed_pca_mgb <- typed_pca_genos %>%
    anti_join(sample_annotation, by = c("id2" = "sample_id")) %>%
    mutate(population = "MGB_biobank") %>%
    unite("sample_id", c("id1", "id2"), sep = "_") %>%
    select(sample_id, population, starts_with("PC"))


typed_pca_df <- bind_rows(typed_pca_mgb, typed_pca_kgp) %>%
    mutate(population = factor(population, levels = names(all_cols))) %>%
    mutate(dataset = ifelse(grepl("MGB", population), "MGB", "1000 Genomes"),
           dataset = factor(dataset, levels = rev(c("1000 Genomes", "MGB")))) 

typed_pca_thresholds <- typed_pca_kgp %>%
    filter(population %in% c("GBR", "CEU", "TSI", "IBS", "FIN")) %>%
    pivot_longer(PC1:PC4, names_to = "pc") %>%
    group_by(pc) %>%
    slice(-which.max(value)) %>%
    slice(-which.min(value)) %>%
    summarise(value = range(value)) %>%
    mutate(i = c("min", "max"),
           value = value * 1.1) %>%
    ungroup() %>%
    pivot_wider(names_from = pc, values_from = value) %>%
    mutate(dataset = "MGB",
           dataset = factor(dataset, levels = c("1000 Genomes", "MGB")))


typed_pca_plot_1 <- typed_pca_df %>%
    ggplot(aes(PC1, PC2, color = population, alpha = dataset)) +
    geom_point(size = .5) +
    scale_color_manual(values = all_cols) +
    scale_alpha_manual(values = c("MGB" = .1, "1000 Genomes" = 1)) +
    facet_wrap(~dataset) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.text = element_text(face = "bold", size = 11),
          panel.spacing = unit(2, "lines")) +
    guides(color = guide_legend(override.aes = list(size = 2), ncol = 2),
           alpha = FALSE)

typed_pca_plot_2 <- typed_pca_df %>%
    ggplot(aes(PC3, PC4, color = population, alpha = dataset)) +
    geom_point(size = .5) +
    scale_color_manual(values = all_cols) +
    scale_alpha_manual(values = c("MGB" = .1, "1000 Genomes" = 1)) +
    facet_wrap(~dataset) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.text = element_text(face = "bold", size = 11),
          panel.spacing = unit(2, "lines")) +
    guides(color = guide_legend(override.aes = list(size = 2), ncol = 2),
           alpha = FALSE)

typed_pca_grid <- plot_grid(typed_pca_plot_1 + theme(legend.position = "none"),
                      NULL,
                      typed_pca_plot_2 + theme(legend.position = "none"),
                      ncol = 1, rel_heights = c(1, .05, 1))

typed_pca_plot <- plot_grid(typed_pca_grid, get_legend(typed_pca_plot_1),
                      nrow = 1, rel_widths = c(1, .3))

ggsave("./plots/pca_genotyped.png", typed_pca_plot, width = 8, height = 5)

typed_mgb_eur_inds <- typed_pca_df %>%
    filter(dataset == "MGB") %>%
    filter(between(PC1, typed_pca_thresholds$PC1[1], typed_pca_thresholds$PC1[2]),
           between(PC2, typed_pca_thresholds$PC2[1], typed_pca_thresholds$PC2[2]),
           between(PC3, typed_pca_thresholds$PC3[1], typed_pca_thresholds$PC3[2]),
           between(PC4, typed_pca_thresholds$PC4[1], typed_pca_thresholds$PC4[2])) %>%
    select(sample_id)

typed_pca_eur_df <- typed_pca_df %>%
    select(sample_id, dataset, population, PC1:PC4) %>%
    mutate(population = as.character(population),
           population = ifelse(sample_id %in% typed_mgb_eur_inds$sample_id, "MGB_eur", population),
           population = factor(population, levels = names(all_cols)))

all_cols[1:2] <- c("grey", "black")

typed_pca_eur_plot_1 <- typed_pca_eur_df %>%
    ggplot(aes(PC1, PC2, color = population, alpha = dataset)) +
    geom_point(size = .5) +
    geom_vline(data = typed_pca_thresholds, 
               aes(xintercept = PC1), linetype = 2, alpha = .5, color = "red") +
    geom_hline(data = typed_pca_thresholds, 
               aes(yintercept = PC2), linetype = 2, alpha = .5, color = "red") +
    scale_color_manual(values = all_cols) +
    scale_alpha_manual(values = c("MGB" = .1, "1000 Genomes" = 1)) +
    facet_wrap(~dataset) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.text = element_text(face = "bold", size = 11),
          panel.spacing = unit(2, "lines")) +
    guides(color = guide_legend(override.aes = list(size = 2), ncol = 2),
           alpha = FALSE)

typed_pca_eur_plot_2 <- typed_pca_eur_df %>%
    ggplot(aes(PC3, PC4, color = population, alpha = dataset)) +
    geom_point(size = .5) +
    geom_vline(data = typed_pca_thresholds, 
               aes(xintercept = PC3), linetype = 2, alpha = .5, color = "red") +
    geom_hline(data = typed_pca_thresholds, 
               aes(yintercept = PC4), linetype = 2, alpha = .5, color = "red") +
    scale_color_manual(values = all_cols) +
    scale_alpha_manual(values = c("MGB" = .1, "1000 Genomes" = 1)) +
    facet_wrap(~dataset) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.text = element_text(face = "bold", size = 11),
          panel.spacing = unit(2, "lines")) +
    guides(color = guide_legend(override.aes = list(size = 2), ncol = 2),
           alpha = FALSE)


typed_pca_eur_grid <- plot_grid(typed_pca_eur_plot_1 + theme(legend.position = "none"),
                          NULL,
                          typed_pca_eur_plot_2 + theme(legend.position = "none"),
                          ncol = 1, rel_heights = c(1, .05, 1))

typed_pca_eur_plot <- plot_grid(typed_pca_eur_grid, get_legend(typed_pca_eur_plot_1),
                          nrow = 1, rel_widths = c(1, .3))

ggsave("./plots/pca_eur_genotyedvars.png", typed_pca_eur_plot, width = 8, height = 5)


# Explore differences among PCAs

genot_eur_inds <- 
    anti_join(typed_pca_eur_df %>%
              filter(population == "MGB_eur") %>%
              select(sample_id),
              pca_eur_df %>%
              filter(population == "MGB_eur") %>%
              select(sample_id)) %>%
    pull(sample_id)

allsnp_eur_inds <- 
    anti_join(pca_eur_df %>%
                  filter(population == "MGB_eur") %>%
                  select(sample_id),
              typed_pca_eur_df %>%
                  filter(population == "MGB_eur") %>%
                  select(sample_id)) %>%
    pull(sample_id)


pca_eur_v2 <- typed_pca_eur_df %>%
    filter(dataset == "MGB") %>%
    mutate(population = case_when(sample_id %in% genot_eur_inds ~ "EUR in genotyped SNPs",
                                  sample_id %in% allsnp_eur_inds ~ "EUR in all SNPs",
                                  TRUE ~ "MGB"))

pc12 <- pca_eur_v2 %>%
    filter(population == "MGB") %>%
    ggplot(aes(PC1, PC2)) +
    geom_point(size = .5, color = "grey") +
    geom_point(data = filter(pca_eur_v2, population != "MGB"),
               aes(PC1, PC2, color = population),
               size = .5) +
    scale_color_manual(values = c("black", "green")) +
    geom_vline(data = typed_pca_thresholds, 
               aes(xintercept = PC1), linetype = 2, alpha = .5, color = "red") +
    geom_hline(data = typed_pca_thresholds, 
               aes(yintercept = PC2), linetype = 2, alpha = .5, color = "red") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(color = NULL) +
    guides(color = guide_legend(override.aes = list(size = 2), ncol = 2))

pc34 <- pca_eur_v2 %>%
    filter(population == "MGB") %>%
    ggplot(aes(PC3, PC4)) +
    geom_point(size = .5, color = "grey") +
    geom_point(data = filter(pca_eur_v2, population != "MGB"),
               aes(PC3, PC4, color = population),
               size = .5) +
    scale_color_manual(values = c("black", "green")) +
    geom_vline(data = typed_pca_thresholds, 
               aes(xintercept = PC3), linetype = 2, alpha = .5, color = "red") +
    geom_hline(data = typed_pca_thresholds, 
               aes(yintercept = PC4), linetype = 2, alpha = .5, color = "red") +
    theme_bw() +
    theme(panel.grid = element_blank())

plot_grid(get_legend(pc12), NULL,
          pc12 + theme(legend.position = "none"), 
          pc34 + theme(legend.position = "none"),
          nrow = 2, 
          rel_heights = c(.1, 1))

ggsave("./plots/pca_difference.png", height = 4)


# LD between Langefeld and Bentham

ld_df <- read_tsv("./sle_variants/sle_ld.tsv") %>%
    mutate(region = paste(chr, bentham),
           region = fct_inorder(region))


ld_plot <- ggplot(ld_df, aes(pos, r2)) +
    geom_hline(yintercept = .6, color = "grey", alpha = .5, linetype = 2) +
    geom_point(aes(color = r2), show.legend = FALSE) +
    geom_text_repel(data = filter(ld_df, r2 >= .6),
                    aes(label = langefeld),
                    size = 3, segment.size = .5, segment.color = "grey45") +
    scale_x_continuous(labels = function(x) ceiling(x/1e6L),
                       breaks = scales::pretty_breaks(3)) +
    scale_y_continuous(breaks = c(0, .5, 1)) +
    scale_color_continuous_tableau() +
    facet_wrap(~region, scales = "free_x") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(hjust = 1, vjust = 1)) +
    coord_cartesian(ylim = c(0, 1.1)) +
    labs(x = "Pos (Mb)", y = expression("r"^2))

ggsave("./plots/sle_ld.png", ld_plot, height = 6)

# Heterozygosity scores
batch_df <- sprintf("./results/eur_females_%s.txt", batches) %>%
    setNames(batches) %>%
    map_df(~tibble(sample_id = read_lines(.)), .id = "batch")

scores_df <- read_tsv("./sle_variants/scores.tsv") %>%
    mutate(subject_id = str_extract(sample_id, "\\d+$")) %>%
    left_join(batch_df, by = "sample_id") %>%
    filter(batch != "0410")
    
candidate_inds <- scores_df %>%
    arrange(desc(het_score), desc(het_score_wt)) %>%
    slice(1:2000)

candidate_inds %>%
    pull(subject_id) %>%
    as.numeric() %>%
    sort() %>%
    write_lines("./candidate_individuals_ids.txt")

scores_plot_df <- scores_df %>%
    mutate(candidate = sample_id %in% candidate_inds$sample_id)

scores_plot_df %>%
    group_by(het_score) %>%
    mutate(y = 1:n()) %>%
    arrange(het_score) %>%
    ggplot(aes(het_score, y = y, color = candidate)) +
    geom_jitter(size = .1, alpha = .25, show.legend = FALSE) +
    scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "purple4")) +
    theme_minimal() +
    theme(panel.grid = element_blank()) +
    labs(x = "Number of heterozygous SLE variants",
         y = "N")

ggsave("./plots/het_scores_dist.png", width = 5, height = 3)

# Duplicates

dups_df <- read_tsv("./sle_variants/scores.tsv") %>%
    mutate(subject_id = str_extract(sample_id, "\\d+$")) %>%
    left_join(batch_df, by = "sample_id") %>%
    add_count(subject_id) %>%
    filter(n > 1) %>%
    select(sample_id, subject_id, het_score) %>%
    left_join(select(het, batch, sample_id = INDV)) %>%
    arrange(subject_id, batch) %>%
    group_by(subject_id) %>%
    mutate(i = 1:n()) %>%
    ungroup() %>%
    select(-sample_id, -batch) %>%
    pivot_wider(names_from = i, values_from = het_score) %>%
    rename(`previous batch` = `1`, `batch 0410` = `2`)

dups_plot <- ggplot(dups_df, aes(`batch 0410`, `previous batch`)) +
    geom_point(size = .25) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          text = element_text(size = 8),
          plot.title = element_text(size = 6)) +
    coord_fixed() +
    labs(title = "Number of SLE variants at which an individual is heterozygote\nfor duplicates in MGB Biobank")

ggsave("./plots/duplicates.png", dups_plot, width = 3, height = 3)

# Heterozygosity per batch

ggplot(scores_df, aes(batch, het_score)) +
    geom_violin() +
    geom_quasirandom(method = "smiley", alpha = .2, size = .1) +
    theme_bw()



## Legacy code #################################################################

ggplot(scores_df, aes(het_score, het_score_wt)) +
    geom_point(aes(color = candidate), size = .5, alpha = .25) +
    scale_color_manual(values = c("black", "green4")) +
    theme_bw()  +
    theme(legend.position = c(.9, .2)) +
    labs(x = "Heterozigosity score", 
         y = "Weighted heterozygosity score",
         color = "Candidate for\nselection?")


ggsave("./plots/het_scores_points.png", width = 6, height = 4)


score_colors <- 
    bind_rows(unweighted = top_n(scores_df, 12, het_score), 
              weighted = top_n(scores_df, 12, het_score_wt), .id = "score") %>%
    select(-top_eur) %>%
    group_by(sample_id) %>%
    mutate(score = ifelse(n_distinct(score) > 1, "both", score))


scores_df %>%
    pivot_longer(het_score:het_score_wt, names_to = "score") %>%
    mutate(score = recode(score, "het_score" = "Unweighted",
                          "het_score_wt" = "Weighted")) %>%
    ggplot(aes(value, fill = top_eur)) +
    geom_density(alpha = .5, size = .1) +
    scale_fill_manual(values = c("TRUE" = "black", "FALSE" = "grey70"),
                      labels = c("Low EUR ancestry", "High EUR ancestry")) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = "top",
          text = element_text(size = 8)) +
    facet_wrap(~score, scales = "free") +
    labs(fill = NULL, x = "Heterozygosity score")

ggsave("./plots/het_score_density.png", width = 4.5, height = 2.5)



## color PCA by heterozygosity

pca_score_df <- pca_df %>%
    filter(dataset == "MGB") %>%
    left_join(scores_df) %>%
    select(sample_id, PC1, PC2, PC3, PC4, het_score, het_score_wt)

pca_score_df %>%
    ggplot(aes(PC1, PC2, color = het_score_wt, fill = het_score_wt)) +
    geom_point(size = .1) +
    scale_color_viridis_c(option = "inferno") +
    scale_fill_viridis_c(option = "inferno", 
                         guide = guide_colorbar(barwidth = 0.5)) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    guides(color = FALSE) +
    labs(fill = "Weighted\nscore")

ggsave("./plots/pca_het_scores.png", width = 4, height = 2.5)





# kmeans

pca_matrix <- pca_df %>%
    select(PC1:PC4) %>%
    as.matrix()

add_cluster <- function(k) {
    
    clusters <- kmeans(pca_matrix, k)$cluster
    
    pca_df %>%
        add_column(cluster = clusters)
}

clusters_df <- map_df(1:10, add_cluster, .id = "k") %>%
    mutate(k = as.integer(k)) %>%
    arrange(k)

clusters_continent <- clusters_df %>%
    filter(population != "MGB_biobank") %>%
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
    labs(x = "K", y = "Homogeneity")

ggsave("./plots/kmeans_accuracy.png", height = 3, width = 5)

clusters_continent %>%
    filter(k == 7, continent == "EUR")

### Illustration (slides)

k5 <- kmeans(pca_matrix, 5)

k5_df <- pca_df %>%
    add_column(cluster = k5$cluster)

k5_cts <- as.data.frame(k5$centers)

pa <- ggplot(k5_df, aes(PC1, PC2)) +
    geom_point(alpha = .05) +
    theme_bw()

pb <- ggplot(k5_df, aes(PC1, PC2)) +
    geom_point(alpha = .05) +
    geom_point(data = k5_cts, 
               aes(PC1, PC2, color = factor(1:5)), 
               size = 4) +
    scale_color_nejm(guide = NULL) +
    theme_bw() +
    theme(panel.grid = element_blank())

pc <- ggplot(k5_df, aes(PC1, PC2, color = factor(cluster))) +
    geom_point(alpha = .2) +
    scale_color_nejm(guide = NULL) +
    theme_bw() +
    theme(panel.grid = element_blank())

plot_grid(pa, pb, pc, ncol = 1)
ggsave("./plots/kmeans_illustration.png", width = 4, height = 6)

###

# Select MGB individuals with high EUR ancestry
continent_colors <- all_cols[c("YRI", "TSI", "PJL", "CDX", "PEL")] %>%
    setNames(levels(clusters_continent$continent)) %>%
    c("MGB" = "Black")

clusters_plot_df <- clusters_df %>%
    left_join(select(clusters_continent, -p), by = c("k", "cluster")) %>%
    filter(dataset == "MGB") %>%
    mutate(continent = as.character(continent),
           continent = ifelse(is.na(continent), "MGB", continent),
           continent = factor(continent, levels = names(continent_colors)))

cluster_plot_1 <- ggplot(clusters_plot_df, aes(PC1, PC2, color = continent)) +
    geom_point(size = .25, alpha = .25) +
    scale_color_manual(values = continent_colors) +
    facet_wrap(~k, nrow = 1) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          strip.text = element_text(face = "bold", size = 11),
          panel.spacing = unit(0, "lines"),
          axis.text = element_blank(),
          legend.position = "top") +
    guides(color = guide_legend(nrow = 1, override.aes = list(size = 2, alpha = 1)))

cluster_plot_2 <- ggplot(clusters_plot_df, aes(PC3, PC4, color = continent)) +
    geom_point(size = .25, alpha = .25) +
    scale_color_manual(values = continent_colors) +
    facet_wrap(~k, nrow = 1) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          strip.text = element_text(face = "bold", size = 11),
          panel.spacing = unit(0, "lines"),
          axis.text = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))

cluster_plot_grid <- 
    plot_grid(cluster_plot_1 + theme(legend.position = "none"),
              cluster_plot_2 + theme(legend.position = "none"),
              ncol = 1)

cluster_plot_final <-
    plot_grid(get_legend(cluster_plot_1), cluster_plot_grid, 
              ncol = 1, rel_heights = c(.1, 1))

ggsave("./plots/clusters.png", cluster_plot_final, height = 4, width = 8)


mgb_eur <- clusters_df %>%
    filter(dataset == "MGB", k == 8, cluster %in% c(4L, 8L)) %>%
    pull(sample_id)




# OR

sle_vars <- read_tsv("./mgb_biobank/sle_variants/sle.tsv")

ggplot(sle_vars, aes(or, log(or))) +
    geom_point()

# Admixture

## Cross-validations to choose best value for K
k_df <- 
    "grep -h CV ./mgb_biobank/results/allchr.merged.pruned.1000G.cv*.log" %>%
    system(intern = TRUE) %>%
    str_remove("CV error ") %>%
    tibble(tmp = .) %>%
    separate(tmp, c("K", "error"), sep = " ", convert = TRUE) %>%
    mutate(K = parse_number(K))

ggplot(k_df, aes(K, error)) +
    geom_point() +
    geom_line(aes(group = 1)) +
    geom_label_repel(aes(label = error), size = 2.5,
                     hjust = "inward", vjust = "inward")  +
    scale_x_continuous(breaks = 1:10) +
    scale_y_continuous(breaks = c(.3, .33, .36)) +
    theme(panel.grid.minor = element_blank()) +
    labs(y = "Cross-validation error")

ggsave("./plots/admixture_cv.png", height = 3)

## Ancestry proportions
ids_1000g <- "./mgb_biobank/results/allchr.merged.pruned.1000G.fam" %>%
    read_delim(delim = " ", col_names = FALSE) %>%
    pull(1)

ancestry_1000g <- 
    "./mgb_biobank/results/allchr.merged.pruned.1000G.%d.Q" %>%
    sprintf(1:10) %>%
    setNames(1:10) %>%
    map_df(~read_delim(., delim = " ", col_names = FALSE) %>%
               add_column(id = ids_1000g, .before = 1) %>%
               pivot_longer(-id, names_to = "cluster", values_to = "q"), 
           .id = "K") %>%
    left_join(sample_annotation, c("id" = "sample_id")) %>%
    mutate(continent = case_when(population %in% names(afr_cols) ~ "AFR",
                                 population %in% names(eur_cols) ~ "EUR",
                                 population %in% names(sas_cols) ~ "SAS",
                                 population %in% names(eas_cols) ~ "EAS",
                                 population %in% names(amr_cols) ~ "AMR")) %>%
    mutate(continent = factor(continent,
                              levels = c("AFR", "EUR", "SAS", "EAS", "AMR")),
           K = as.integer(K),
           cluster = factor(cluster, levels = sprintf("X%d", 1:10))) %>%
    arrange(continent, population, cluster, id) %>%
    mutate(population = fct_inorder(population),
           id = fct_inorder(id))


admix1000G <- ggplot(ancestry_1000g, aes(q, id, fill = cluster)) +
    geom_bar(stat = "identity", position = "fill", width = 1) +
    facet_grid(continent~K, scales = "free") +
    scale_fill_npg() +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.spacing = unit(0, "lines"),
          legend.position = "none") +
    labs(y = NULL)

ggsave("./plots/admixture.1000G.allK.png", admix1000G, height = 7)


## MGB biobank ancestry
ids_1000g <- "./mgb_biobank/results/allchr.refpanel.fam" %>%
    read_delim(delim = " ", col_names = FALSE) %>%
    pull(X1)

cluster_identity <- "./mgb_biobank/results/allchr.refpanel.3.Q" %>%
    read_delim(delim = " ", col_names = FALSE) %>%
    add_column(id = ids_1000g, .before = 1) %>%
    left_join(sample_annotation, by = c("id" = "sample_id")) %>%
    pivot_longer(X1:X3, names_to = "cluster", values_to = "q") %>%
    group_by(population, cluster) %>%
    summarise(q = mean(q)) %>%
    slice(which.max(q)) %>%
    ungroup() %>%
    select(cluster, population)


ids_mgb <- "./mgb_biobank/results/allchr.merged.pruned.MGB.fam" %>%
    read_delim(delim = " ", col_names = FALSE) %>%
    pull(2)

ancestry_mgb <- "./mgb_biobank/results/allchr.merged.pruned.MGB.3.Q" %>%
    read_delim(delim = " ", col_names = FALSE) %>%
    add_column(id = ids_mgb, .before = 1) %>%
    pivot_longer(X1:X3, names_to = "cluster", values_to = "q") %>%
    group_by(id) %>%
    mutate(likely_group = cluster[which.max(q)],
           group_q = max(q)) %>%
    ungroup() %>%
    arrange(likely_group, desc(group_q)) %>%
    mutate(id = fct_inorder(id))

refpanel_colors <- 
    #all_cols[c("CDX", "YRI", "GBR", "GIH")] %>%
    all_cols[c("TSI", "CHS", "YRI")] %>%
    setNames(paste0("X", 1:3))

ggplot(ancestry_mgb, aes(id, q, fill = cluster), color = NULL) +
    geom_col(width = 1) +
    scale_fill_manual(values = refpanel_colors) +
    scale_y_continuous(breaks = c(0, .5, 1), labels = scales::percent) +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid = element_blank(), 
          strip.text = element_blank(),
          panel.spacing = unit(.25, "lines"),
          legend.position = "none") +
    facet_grid(~likely_group, scales = "free", space = "free") +
    labs(x = "Individual", y = "Ancestry %")

ggsave("./plots/admixture_mgb_k3.png")

