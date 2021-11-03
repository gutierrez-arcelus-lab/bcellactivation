library(tidyverse)
library(ggsci)
library(ggrepel)


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


#PCA

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

