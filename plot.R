library(tidyverse)
library(ggsci)
library(ggrepel)
library(RColorBrewer)
library(cowplot)

# Most expressed biotypes
gene_df <- read_tsv("./data/gene_quants.tsv")

transc_df <- read_tsv("./data/transcript_quants.tsv")

transc_types_summary <- transc_df %>%
    group_by(condition_id, transcript_type) %>%
    summarise(cpm = sum(cpm),
              tpm = sum(tpm)) %>%
    group_by(condition_id) %>%
    mutate(cpm = cpm/sum(cpm),
           tpm = tpm/sum(tpm)) %>%
    group_by(transcript_type) %>%
    mutate(transcript_type = ifelse(any(tpm>0.01), transcript_type, "other")) %>%
    group_by(condition_id, transcript_type) %>%
    summarise(avg_cpm = sum(cpm),
              avg_tpm = sum(tpm)) %>%
    ungroup()

biotypes_order <-
    transc_types_summary %>%
    filter(condition_id == "16hr_resting") %>%
    arrange(-avg_cpm) %>%
    pull(transcript_type) %>%
    fct_inorder() %>%
    fct_relevel("other", after = Inf) %>%
    levels() %>%
    rev()

transc_types_summary_reordered <- transc_types_summary %>%
    mutate(transcript_type = factor(transcript_type, levels = biotypes_order),
           condition_id = factor(condition_id, 
                                 levels = c("16hr_resting", 
                                            "24hr_IgG", "72hr_IgG",
                                            "24hr_RSQ", "72hr_RSQ")))

transc_types_summary_reordered %>%
    select(transcript_type, condition_id, TPM = avg_tpm) %>%
    ggplot(aes(x = condition_id, y = TPM, fill = transcript_type)) +
    geom_col() +
    scale_x_discrete(labels = function(x) sub("_", "\n", x)) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_npg() +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(vjust = 5)) +
    labs(y = "% of total expression", x = NULL, 
         fill = "Transcript biotype")
    
ggsave("./plots/total_expression.png", width = 6)


# Fold changes
logfc <- read_tsv("./data/logfc.tsv") %>%
    mutate(gene_type = ifelse(gene_type %in% biotypes_order, gene_type, "other"),
           gene_type = factor(gene_type, levels = biotypes_order))

ggplot(logfc, aes(cpm, log2fc)) +
    geom_point(alpha = .1, size = .5) +
    facet_wrap(~condition_id, scales = "free", ncol = 2) +
    theme_bw() +
    theme(panel.grid = element_blank())  +
    labs(x = "Counts per million in condition",
         y = expression(paste("Log"[2], FC)))

ggsave("./plots/foldchange.png", width = 5)

logfc %>%
    filter(cpm > 10) %>%
    ggplot(aes(cpm, log2fc)) +
    geom_point(alpha = .1, size = .5) +
    facet_wrap(~condition_id, scales = "free", ncol = 2) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(x = "Counts per million in condition",
     y = expression(paste("Log"[2], FC)),
     title = "CPM > 10")

ggsave("./plots/foldchange_subset.png", width = 5)

logfc %>%
    group_by(condition_id) %>%
    summarise(downregulated = mean(log2fc < 0),
              upregulated = mean(log2fc > 0)) %>%
    pivot_longer(downregulated:upregulated, names_to = "direction") %>%
    ggplot(aes(condition_id, value, fill = direction)) +
    geom_col(position = "dodge") +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_npg() +
    labs(x = NULL, y = "% of genes", fill = NULL)

ggsave("./plots/foldchange_summary.png", width = 5)



#PCA on expression data

pca <- "./data/pca/pheno_pcs.pca" %>%
    read_delim(delim = " ") %>%
    pivot_longer(-1, names_to = "condition_id") %>%
    mutate(pc = str_extract(SampleID, "(PC\\d$)")) %>%
    select(pc, condition_id, value) %>%
    pivot_wider(names_from = pc, values_from = value) %>%
    mutate(condition_id = factor(condition_id, 
                                 levels = c("16hr_resting", "24hr_IgG",
                                            "72hr_IgG", "24hr_RSQ", "72hr_RSQ")))

ggplot(pca, aes(PC1, PC2, color = condition_id)) +
    geom_point(size = 5) +
    scale_color_manual(values = c("16hr_resting" = "grey35",
                                  "24hr_IgG" = "slateblue",
                                  "72hr_IgG" = "cornflowerblue",
                                  "24hr_RSQ" = "tomato3",
                                  "72hr_RSQ" = "salmon")) +
    theme_bw()

ggsave("./plots/pca_bcell_expression.png")

# Correlations between conditions
gene_names <- read_tsv("./data/transc_to_gene.tsv") %>%
    distinct(gene_id, gene_name)

gene_v2_df <- gene_df %>%
    select(-tpm) %>%
    pivot_wider(names_from = condition_id, values_from = cpm) %>%
    pivot_longer(-(1:2), names_to = "condition", values_to = "cpm")

cor_df <- gene_v2_df %>%
    group_by(condition) %>%
    summarise(y = max(cpm),
              x = min(`16hr_resting`),
              rho = cor(cpm, `16hr_resting`, method = "spearman"),
              rho = round(rho, 2),
              rho_label = paste("rho == ", rho)) %>%
    ungroup()

gene_labels_df <- gene_v2_df %>%
    filter(`16hr_resting` > 2500 , cpm > 2500) %>%
    mutate(fc = pmax(cpm, `16hr_resting`)/pmin(cpm, `16hr_resting`),
           fc = replace_na(fc, 0)) %>%
    group_by(condition) %>%
    top_n(10, fc) %>%
    ungroup() %>%
    left_join(gene_names)    

ggplot(gene_v2_df, aes(`16hr_resting`, cpm)) +
    geom_abline() +
    geom_point(alpha = .25) +
    geom_text(data = cor_df, aes(x, y * 1.1, label = rho_label),
              hjust = "inward", vjust = "inward",
              parse = TRUE, size = 3) +
    geom_text_repel(data = gene_labels_df, 
                    aes(label = gene_name),
                    size = 3) +
    facet_wrap(~condition, scales = "free", ncol = 2) +
    theme_bw()

ggsave("./plots/scatter_resting_conditions.png")


condition_df <- logfc %>%
    select(-cpm) %>%
    pivot_wider(names_from = condition_id, values_from = log2fc)

igG_diff <- condition_df %>%
    filter(`24hr_IgG` < -50 & `72hr_IgG` > -20 | 
           `72hr_IgG` < -50 & `24hr_IgG` > -20) %>%
    left_join(gene_names)

rsq_diff <- condition_df %>%
    filter(`24hr_RSQ` < -50 & `72hr_RSQ` > -20 | 
               `72hr_RSQ` < -50 & `24hr_RSQ` > -20) %>%
    left_join(gene_names)

igg <- ggplot(filter(condition_df, ! gene_id %in% igG_diff$gene_id),
       aes(`24hr_IgG`, `72hr_IgG`)) +
    geom_abline(color = "grey", linetype = 2) +
    geom_point(size = .5) +
    geom_point(data = igG_diff, 
               aes(`24hr_IgG`, `72hr_IgG`, color = gene_type),
               size = .5) +
    scale_color_npg() +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(x = expression(paste("Log"[2], "FC (24hr)")),
         y = expression(paste("Log"[2], "FC (72hr)")),
         title = "IgG")

rsq <- ggplot(filter(condition_df, ! gene_id %in% rsq_diff$gene_id),
       aes(`24hr_RSQ`, `72hr_RSQ`)) +
    geom_abline(color = "grey", linetype = 2) +
    geom_point(size = .5) +
    geom_point(data = rsq_diff, 
               aes(`24hr_RSQ`, `72hr_RSQ`, color = gene_type),
               size = .5) +
    scale_color_npg() +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(x = expression(paste("Log"[2], "FC (24hr)")),
         y = expression(paste("Log"[2], "FC (72hr)")),
         title = "RSQ")

plot_grid(igg, rsq, ncol = 1)
ggsave("./plots/fc_24vs72.png", width = 5, height = 4)

# MGB biobank


## PCA on genotype data

### Metadata for 1000G samples
index_1000G <- "https://raw.githubusercontent.com/igsr/1000Genomes_data_indexes/master/data_collections/1000_genomes_project/1000genomes.sequence.index"

sample_annotation <- read_tsv(index_1000G, comment = "##") %>%
    select(sample_id = SAMPLE_NAME,
           population = POPULATION) %>%
    distinct()


### Colors
mgb_cols <- c("black", "cornflowerblue") %>%
    setNames(c("MGB_biobank", "MGB_most_EUR"))

afr_cols <- brewer.pal("Oranges", n = 9)[3:9] %>%
    setNames(c("ACB", "ESN", "GWD", "LWK", "MSL", "YRI", "ASW"))

eur_cols <- brewer.pal("Blues", n =9)[5:9] %>%
    setNames(c("GBR", "IBS", "TSI", "CEU", "FIN"))

sas_cols <- brewer.pal("Greens", n =9)[5:9] %>%
    setNames(c("BEB", "PJL", "GIH", "ITU", "STU"))

eas_cols <- brewer.pal("Purples", n =9)[5:9] %>%
    setNames(c("CHB", "CHS", "CDX", "KHV", "JPT"))

amr_cols <- c("lightpink1", "hotpink", "hotpink3", "deeppink") %>%
    setNames(c("MXL", "CLM", "PEL", "PUR"))

all_cols <- c(mgb_cols, afr_cols, eur_cols, sas_cols, eas_cols, amr_cols)

### read PCA results into R
pca_genos <-
    file.path("/lab-share/IM-Gutierrez-e2/Public/vitor/ase/mgb_biobank/results",
    "allchr.merged.pruned.pca.eigenvec") %>%
    read_table(col_names = FALSE) %>%
    select(-1) %>%
    select(X2:X6) %>%
    setNames(c("sample_id", paste0("PC", 1:4)))

pca_kgp <- pca_genos %>%
    inner_join(sample_annotation) 

pca_mgb <- pca_genos %>%
    anti_join(sample_annotation) %>%
    mutate(population = "MGB_biobank")


pca_df <- bind_rows(pca_mgb, pca_kgp) %>%
    select(sample_id, population, PC1:PC4) %>%
    mutate(population = factor(population, levels = names(all_cols))) %>%
    mutate(dataset = ifelse(grepl("MGB", population), "MGB", "1000 Genomes"),
           dataset = factor(dataset, levels = rev(c("1000 Genomes", "MGB")))) 


pca_plot_1 <- pca_df %>%
    ggplot(aes(PC1, PC2, color = population, alpha = dataset)) +
    geom_point(size = .75) +
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
    geom_point(size = .75) +
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
    select(k, cluster, continent, p)
    
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

# Select MGB individuals with high EUR ancestry
mgb_eur <- clusters_df %>%
    filter(dataset == "MGB", k == 7, cluster %in% c(3L, 7L)) %>%
    pull(sample_id)

continent_colors <- all_cols[c("YRI", "GBR", "BEB", "CHS", "PEL")] %>%
    setNames(levels(clusters_continent$continent))

clusters_plot_df <- clusters_df %>%
    left_join(clusters_continent) %>%
    filter(dataset != "MGB") 

cluster_plot_1 <- ggplot(clusters_plot_df, aes(PC1, PC2, color = continent)) +
    geom_point(size = .75, alpha = .5) +
    scale_color_manual(values = continent_colors) +
    facet_wrap(~k, nrow = 1) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          strip.text = element_text(face = "bold", size = 11),
          panel.spacing = unit(0, "lines"),
          axis.text = element_blank(),
          legend.position = "top") +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))

cluster_plot_2 <- ggplot(clusters_plot_df, aes(PC3, PC4, color = continent)) +
    geom_point(size = .75, alpha = .5) +
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


# Heterozygosity scores

scores_df <- read_tsv("./mgb_biobank/sle_variants/scores.tsv") %>%
    extract(sample_id, c("sample_id"), ".+_(.+)") %>%
    mutate(top_eur = sample_id %in% mgb_eur)

scores_eur_df <- scores_df %>%
    filter(top_eur == TRUE)


ggplot(scores_df, aes(het_score, het_score_wt)) +
    geom_jitter(size = .75, alpha = .1) +
    theme_bw()  +
    labs(x = "Heterozigosity score", y = "Weighted heterozygosity score")

ggsave("./plots/het_scores_jitterplot.png", width = 5, height = 4)


scores_df %>%
    pivot_longer(het_score:het_score_wt, names_to = "score") %>%
    mutate(score = recode(score, "het_score" = "unweighted",
                          "het_score_wt" = "weighted")) %>%
    ggplot(aes(value, color = top_eur, fill = top_eur)) +
    geom_density(alpha = .5, size = .5) +
    theme_bw() +
    scale_color_manual(values = c("TRUE" = "midnightblue",
                                  "FALSE" = "grey70"), guide = "none") +
    scale_fill_manual(values = c("TRUE" = "midnightblue",
                                 "FALSE" = "grey70")) +
    facet_wrap(~score, scales = "free") +
    labs(fill = "Top\nEuropean\nAncestry?")
    
ggsave("./plots/het_score_density.png", width = 6)


# OR

sle_vars <- read_tsv("./mgb_biobank/sle_variants/sle.tsv")

ggplot(sle_vars, aes(or, log(or))) +
    geom_point()

# Admixture

## Cross-validations to choose best value for K
k_df <- system("grep -h CV ./mgb_biobank/results/admix.1000G.cv*.log", intern = TRUE) %>%
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
ids_1000g <- "./mgb_biobank/results/allchr.1000G.fam" %>%
    read_delim(delim = " ", col_names = FALSE) %>%
    pull(1)

ancestry_1000g <- 
    "./mgb_biobank/results/allchr.1000G.%d.Q" %>%
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



## Plot results for reference panel
ids_refpanel <- "./mgb_biobank/results/allchr.refpanel.fam" %>%
    read_delim(delim = " ", col_names = FALSE) %>%
    pull(1)

ancestry_refpanel <- "./mgb_biobank/results/allchr.refpanel.3.Q" %>%
    read_delim(delim = " ", col_names = FALSE) %>%
    add_column(id = ids_refpanel, .before = 1) %>%
    left_join(sample_annotation, c("id" = "sample_id")) %>%
    mutate(continent = case_when(population %in% names(afr_cols) ~ "AFR",
                                 population %in% names(eur_cols) ~ "EUR",
                                 population %in% names(sas_cols) ~ "SAS",
                                 population %in% names(eas_cols) ~ "EAS",
                                 population %in% names(amr_cols) ~ "AMR")) %>%
    mutate(continent = factor(continent, 
                              levels = c("EUR", "AFR", "SAS", "EAS", "AMR"))) %>%
    pivot_longer(X1:X3, names_to = "cluster", values_to = "q") %>%
    arrange(continent, population, id) %>%
    mutate(population = fct_inorder(population))

# refpanel_colors <- c("X1" = all_cols[["TSI"]], "X2" = all_cols[["YRI"]],
#                      "X3" = all_cols[["BEB"]], "X4" = all_cols[["CHS"]],
#                      "X5" = all_cols[["CLM"]])

refpanel_colors <- c("X1" = all_cols[["TSI"]], "X2" = all_cols[["YRI"]],
                     "X3" = all_cols[["CHB"]])

ggplot(ancestry_refpanel, aes(id, q, fill = cluster)) +
    geom_bar(stat = "identity", position = "fill", width = 1) +
    facet_wrap(~population, scales = "free", nrow = 1) +
    scale_fill_manual(values = refpanel_colors) +
    scale_y_continuous(breaks = c(0, .5, 1)) +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid = element_blank()) +
    labs(x = NULL)

ggsave("./plots/admixture_refpanel.png")


## MGB biobank ancestry
ids_mgb <- "./mgb_biobank/results/allchr.MGB.fam" %>%
    read_delim(delim = " ", col_names = FALSE) %>%
    pull(1)

ancestry_mgb <- "./mgb_biobank/results/allchr.MGB.3.Q" %>%
    read_delim(delim = " ", col_names = FALSE) %>%
    add_column(id = ids_mgb, .before = 1) %>%
    arrange(X1) %>%
    mutate(id = fct_inorder(id)) %>%
    pivot_longer(X1:X3, names_to = "cluster", values_to = "q")

ggplot(ancestry_mgb, aes(id, q, fill = cluster)) +
    geom_bar(stat = "identity", position = "fill", width = 1) +
    scale_fill_manual(values = refpanel_colors) +
    scale_y_continuous(breaks = c(0, .5, 1)) +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid = element_blank()) +
    labs(x = NULL)

ggsave("./plots/admixture_mgb_k3.png")
