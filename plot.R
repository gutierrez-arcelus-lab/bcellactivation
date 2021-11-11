library(tidyverse)
library(ggsci)
library(ggrepel)
library(RColorBrewer)


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
mgb_cols <- "grey70" %>%
    setNames("MGB_biobank")

afr_cols <- brewer.pal("Oranges", n = 9)[4:9] %>%
    c("black") %>%
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

pca_mgb <- pca_genos %>%
    anti_join(sample_annotation) %>%
    mutate(population = "MGB_biobank")

pca_df <- bind_rows(pca_mgb, pca_kgp) %>%
    select(sample_id, population, PC1:PC3) %>%
    mutate(population = factor(population, levels = names(all_cols)))

pca_for_plot <- 
    bind_rows(
        select(pca_df, sample_id, population, x = PC1, y = PC2) %>%
            mutate(comparison = "PC1 vs PC2"),
        select(pca_df, sample_id, population, x = PC2, y = PC3) %>%
            mutate(comparison = "PC2 vs PC3"))


sizes <- ifelse(names(all_cols) == "MGB_biobank", 2, 1) %>%
    setNames(names(all_cols))

ggplot(pca_for_plot, aes(x, y, color = population, size = population)) +
    geom_point() +
    scale_color_manual(values = all_cols) +
    scale_size_manual(values = sizes) +
    facet_wrap(~comparison, ncol = 1, scales = "free") +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          strip.text = element_text(face = "bold", size = 11),
          panel.spacing = unit(2, "lines")) +
    guides(color = guide_legend(override.aes = list(size = 2),
                                ncol = 2)) +
    labs(x = NULL, y = NULL)

ggsave("./plots/pca.png", height = 6)

k_df <- system("grep -h CV ./mgb_biobank/results/admix.cv*.log", intern = TRUE) %>%
    str_remove("CV error ") %>%
    tibble(tmp = .) %>%
    separate(tmp, c("K", "error"), sep = " ", convert = TRUE) %>%
    mutate(K = parse_number(K))

ggplot(k_df, aes(K, error)) +
    geom_point() +
    geom_line(aes(group = 1)) +
    labs(y = "Cross-validation error")


