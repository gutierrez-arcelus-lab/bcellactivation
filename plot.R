library(tidyverse)
library(ggsci)


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
    facet_wrap(~condition_id, scales = "free", ncol = 1) +
    theme_bw() +
    theme(panel.grid = element_blank()) 

#PCA

pca <- "./data/pca/pheno_pcs.pca" %>%
    read_delim(delim = " ") %>%
    pivot_longer(-1, names_to = "condition_id") %>%
    mutate(pc = str_extract(SampleID, "(PC\\d$)")) %>%
    select(pc, condition_id, value) %>%
    pivot_wider(names_from = pc, values_from = value)

ggplot(pca, aes(PC1, PC2, color = condition_id)) +
    geom_point(size = 5) +
    scale_color_viridis_d() +
    theme_bw()

# Correlations between conditions

genes_df <- read_tsv("./phenotypes.bed.gz") %>%
    select(gene_id = gid, contains("hr_")) %>%
    pivot_longer(-(1:2), names_to = "condition_id", values_to = "tpm") %>%
    select(gene_id, resting = `16hr_resting`, condition_id, tpm) %>%
    arrange(gene_id, condition_id)





