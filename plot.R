library(tidyverse)
library(ggsci)


# Most expressed biotypes
quants_df <- read_tsv("./transcript_quants.tsv")

biotypes_summary <- quants_df %>%
    group_by(id, transcript_biotype) %>%
    summarise(S = sum(tpm)) %>%
    group_by(id) %>%
    mutate(S = S/sum(S)) %>%
    group_by(transcript_biotype) %>%
    mutate(transcript_biotype = ifelse(any(S>0.01), transcript_biotype, "other")) %>%
    group_by(id, transcript_biotype) %>%
    summarise(S = sum(S)) %>%
    ungroup()

biotypes_order <-
    biotypes_summary %>%
    filter(id == "16hr_resting") %>%
    arrange(-S) %>%
    pull(transcript_biotype) %>%
    fct_inorder() %>%
    fct_relevel("other", after = Inf) %>%
    levels() %>%
    rev()

biotypes_summary_reordered <- biotypes_summary %>%
    mutate(transcript_biotype = factor(transcript_biotype, levels = biotypes_order),
           id = factor(id, levels = c("16hr_resting", "24hr_IgG", "72hr_IgG",
                                      "24hr_RSQ", "72hr_RSQ")))

ggplot(biotypes_summary_reordered, aes(id, S, fill = transcript_biotype)) +
    geom_bar(stat = "identity", position = "fill") +
    labs(x = NULL, y = "Proportion of total expression") +
    scale_fill_npg() +
    theme_minimal() +
    theme(panel.grid = element_blank())

ggsave("./plots/transcript_biotypes.png")


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





