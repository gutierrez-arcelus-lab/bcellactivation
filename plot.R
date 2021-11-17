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
    rename(CPM = avg_cpm, TPM = avg_tpm) %>%
    pivot_longer(CPM:TPM, names_to = "est") %>%
    ggplot(aes(condition_id, value, fill = transcript_type)) +
    geom_bar(stat = "identity", position = "fill") +
    facet_wrap(~est, scales="free") +
    labs(x = NULL, y = "Proportion of total expression") +
    scale_fill_npg() +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 1.25, hjust = 1))

ggsave("./plots/transcript_biotypes.png")


# Fold changes
logfc <- read_tsv("./data/logfc.tsv") %>%
    mutate(gene_type = ifelse(gene_type %in% biotypes_order, gene_type, "other"),
           gene_type = factor(gene_type, levels = biotypes_order))

ggplot(logfc, aes(cpm, log2fc)) +
    geom_point(alpha = .1, size = .5) +
    facet_wrap(~condition_id, scales = "free", ncol = 2) +
    theme_bw() +
    theme(panel.grid = element_blank()) 

ggsave("./plots/fc.png")

ggplot(logfc, aes(cpm, log2fc)) +
    geom_point(alpha = .1, size = .5) +
    facet_wrap(~condition_id, scales = "free", ncol = 2) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    coord_cartesian(xlim = c(0, 200), ylim = c(-6, 6))

ggsave("./plots/fc_subset.png")


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
    "/lab-share/IM-Gutierrez-e2/Public/vitor/ase/mgb_biobank/results/plink_pca.eigenvec" %>%
    read_table2(col_names = FALSE) %>%
    select(-1) %>%
    select(X2:X12) %>%
    setNames(c("sample_id", paste0("PC", 1:10)))

pca_kgp <- pca_genos %>%
    inner_join(sample_annotation) 

mgb_mosteur <- read_lines("./mgb_biobank/mgb_mosteuropean.txt")

pca_mgb <- pca_genos %>%
    anti_join(sample_annotation) %>%
    mutate(population = ifelse(sample_id %in% mgb_mosteur, "MGB_most_EUR", "MGB_biobank"))

pca_df <- bind_rows(pca_mgb, pca_kgp) %>%
    select(sample_id, population, PC1:PC3) %>%
    mutate(population = factor(population, levels = names(all_cols)))

pca_for_plot <- 
    bind_rows(
        select(pca_df, sample_id, population, x = PC1, y = PC2) %>%
            mutate(comparison = "PC1 vs PC2"),
        select(pca_df, sample_id, population, x = PC2, y = PC3) %>%
            mutate(comparison = "PC2 vs PC3")) %>%
    mutate(dataset = ifelse(grepl("MGB", population), "MGB", "1000 Genomes"))

ggplot(pca_for_plot, aes(x, y, color = population, alpha = dataset)) +
    geom_point(size = .75) +
    scale_color_manual(values = all_cols) +
    scale_alpha_manual(values = c("MGB" = .1, "1000 Genomes" = 1)) +
    facet_grid(comparison~dataset, scales = "free") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.text = element_text(face = "bold", size = 11),
          panel.spacing = unit(2, "lines")) +
    guides(color = guide_legend(override.aes = list(size = 2), ncol = 2),
           alpha = FALSE) +
    labs(x = NULL, y = NULL)

ggsave("./plots/pca.png", width = 8, height = 5)

# Heterozygosity scores

scores_df <- read_tsv("./mgb_biobank/sle_variants/scores.tsv") %>%
    extract(sample_id, c("sample_id"), ".+_(.+)")

top_scores <- 
    bind_rows(het = scores_df %>% top_n(12, het_score),
              het_wt = scores_df %>% top_n(12, het_score_wt),
              .id = "score")

ggplot(scores_df, aes(het_score, het_score_wt)) +
    geom_jitter()

top_scores %>%
    select(score, sample_id) %>%
    group_by(sample_id) %>%
    mutate(shared = all(c("het", "het_wt") %in% score)) %>%
    ggplot(aes(score, sample_id, color = shared)) +
    geom_point() +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 4)) +
    scale_color_npg() +
    labs(y = NULL)
    
ggsave("./plots/het_score.png", width = 4)




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

ancestry_1000g <- "./mgb_biobank/results/allchr.1000G.9.Q" %>%
    read_delim(delim = " ", col_names = FALSE) %>%
    add_column(id = ids_1000g, .before = 1) %>%
    left_join(sample_annotation, c("id" = "sample_id")) %>%
    mutate(continent = case_when(population %in% names(afr_cols) ~ "AFR",
                                 population %in% names(eur_cols) ~ "EUR",
                                 population %in% names(sas_cols) ~ "SAS",
                                 population %in% names(eas_cols) ~ "EAS",
                                 population %in% names(amr_cols) ~ "AMR")) %>%
    mutate(continent = factor(continent,
                              levels = c("AFR", "EUR", "SAS", "EAS", "AMR"))) %>%
    pivot_longer(X1:X9, names_to = "cluster", values_to = "q") %>%
    arrange(continent, population, id) %>%
    mutate(population = fct_inorder(population))

ancestry_colors <- c("X1" = amr_cols[["CLM"]], "X2" = afr_cols[["YRI"]],
                     "X3" = eur_cols[["IBS"]], "X4" = sas_cols[["ITU"]],
                     "X5" = eas_cols[["CDX"]])

ggplot(ancestry_1000g, aes(id, q, fill = cluster)) +
    geom_bar(stat = "identity", position = "fill", width = 1) +
    facet_wrap(~population, scales = "free", ncol = 4) +
    #scale_fill_manual(values = ancestry_colors) +
    scale_y_continuous(breaks = c(0, .5, 1)) +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid = element_blank()) +
    labs(x = NULL)

ggsave("./plots/admixture.1000G.png", height = 5)

# 
# ## select individuals for a reduced reference panel
# ancestry_1000g %>%
#     filter(q > .975) %>%
#     filter(population %in% c("YRI", "STU", "CEU", "CHS", "PEL")) %>%
#     pull(id) %>%
#     {paste(., ., sep = "_")} %>%
#     write_lines("./mgb_biobank/ref_panel_ids.txt")


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
