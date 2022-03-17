library(tidyverse)
library(ggsci)
library(ggthemes)
library(ggrepel)
library(RColorBrewer)
library(cowplot)

# Most expressed biotypes
gene_df <- read_tsv("./data/gene_quants.tsv")

gene_names <- read_tsv("./data/transc_to_gene.tsv") %>%
     distinct(gene_id, gene_name)

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

transc_summary <- transc_types_summary %>%
    mutate(transcript_type = factor(transcript_type, levels = biotypes_order),
           condition_id = factor(condition_id, 
                                 levels = c("16hr_resting", 
                                            "24hr_IgG", "72hr_IgG",
                                            "24hr_RSQ", "72hr_RSQ")))

transc_summary %>%
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



# Genes in Bentham et al.
library(rvest)

bentham_genes <- 
    "https://www.nature.com/articles/ng.3434/tables/2" %>%
    read_html() %>%
    html_node("table") %>%
    html_table(header = TRUE, fill = TRUE) %>%
    select(gene = `Likely causal genec`) %>%
    slice(-1) %>%
    separate_rows(gene, sep = ",") %>%
    mutate(gene = trimws(gene)) %>%
    filter(gene != "") %>%
    pull(gene)

bentham_genes[bentham_genes == "CXorf21"] <- "TASL"

conditions <- c("16hr_resting", "24hr_IgG", "72hr_IgG", "24hr_RSQ", "72hr_RSQ")

bentham_tpm <- gene_df %>%
    left_join(gene_names) %>%
    filter(gene_name %in% bentham_genes) %>%
    mutate(condition_id = factor(condition_id, levels = conditions))

bentham_plot <- ggplot(bentham_tpm, aes(gene_name, tpm)) +
    geom_col(aes(fill = condition_id), position = "dodge",
             color = "black", size = .15) +
    scale_fill_manual(values = c("grey70", "gold2", "gold3",
                                 "mediumpurple1", "mediumpurple3"),
                      labels = function(x) sub("_", "\n", x)) +
    scale_y_continuous(breaks = scales::pretty_breaks(3)) +
    facet_wrap(~gene_name, scales = "free") +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_blank(),
          legend.position = c(.75, .05),
          legend.direction = "horizontal",
          legend.background = element_blank(),
          legend.box.background = element_rect(color = "black")) +
    labs(x = NULL, y = "TPM", fill = NULL)

ggsave("./plots/bentham.png", bentham_plot, height = 6, width = 8)


bentham_fc <- bentham_tpm %>%
    mutate(gene_name = as.character(gene_name)) %>%
    select(-cpm) %>%
    mutate(tpm = tpm + 1L) %>%
    pivot_wider(names_from = condition_id, values_from = tpm) %>%
    pivot_longer(-(1:3), names_to = "condition_id", values_to = "tpm") %>%
    select(gene_id, gene_name, condition_id, resting_tpm = 3, tpm) %>%
    mutate(fc = log2(tpm) - log2(resting_tpm)) %>%
    group_by(condition_id, gene_id) %>%
    filter(resting_tpm > 10, tpm > 10) %>%
    filter(fc > 0.5 | fc < -0.5) %>%
    ungroup() %>%
    mutate(condition_id = factor(condition_id, levels = conditions))

library(tidytext)

bentham_fc_plot <- ggplot(bentham_fc, aes(fc, reorder_within(gene_name, fc, condition_id))) +
    geom_vline(xintercept = -3:3, linetype = 2, color = "grey", alpha = .33) +
    geom_col(aes(fill = fc), show.legend = FALSE) +
    facet_wrap(~condition_id, scales = "free_y") +
    scale_y_reordered() +
    scale_fill_gradient2() +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          plot.caption = element_text(hjust = 0),
          axis.text.y = element_text(size = 7, face = "bold")) +
    labs(x = expression(paste("Log"[2], "FC TPM")),
         y = NULL,
         title = "Fold change in respect to resting state\nfor genes in Bentham et al.",
         caption = "Genes with TPM > 10 in both resting and stim and |FC| > 0.5") 
    
ggsave("./plots/bentham_fc.png", bentham_fc_plot, height = 6.5)



# TLRs and IG genes

immune_genes <- read_tsv("./data/transc_to_gene.tsv") %>%
  distinct(gene_id, gene_name, gene_type) %>%
  filter((grepl("TLR", gene_name) & gene_type == "protein_coding") |
           grepl("^IG_.+_gene$", gene_type) |
           gene_name == "UNC93B1")
    

immune_df <- gene_df %>%
  inner_join(immune_genes, by = "gene_id") %>%
  mutate(gene_name = case_when(grepl("^IG", gene_type) ~ gene_type,
                               TRUE ~ gene_name)) %>%
  group_by(condition_id, gene_name) %>%
  summarise(tpm = sum(tpm)) %>%
  ungroup() %>%
  mutate(condition_id = factor(condition_id, levels = conditions)) %>%
  group_by(gene_name) %>%
  filter(any(tpm > 5)) %>%
  ungroup() %>%
  mutate(gene_name = factor(gene_name, levels = immune_gene_levels))

ggplot(immune_df, aes(gene_name, tpm)) +
  geom_col(aes(fill = condition_id), position = "dodge",
           color = "black", size = .15) +
  scale_fill_manual(values = c("grey70", "gold2", "gold3",
                               "mediumpurple1", "mediumpurple3")) +
  scale_y_continuous(breaks = scales::pretty_breaks(3)) +
  facet_wrap(~gene_name, scales = "free", ncol = 3) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank()) +
  labs(x = NULL, y = "TPM", fill = NULL)

ggsave("./plots/immune_genes.png", width = 6, height = 4)


ig_c_genes <- gene_df %>%
  inner_join(immune_genes, by = "gene_id") %>%
  filter(gene_type == "IG_C_gene") %>%
  mutate(condition_id = factor(condition_id, levels = conditions)) %>%
  group_by(gene_name) %>%
  filter(any(tpm > 5)) %>%
  ungroup()

ggplot(ig_c_genes, aes(gene_name, tpm)) +
  geom_col(aes(fill = condition_id), position = "dodge",
           color = "black", size = .15) +
  scale_fill_manual(values = c("grey70", "gold2", "gold3",
                               "mediumpurple1", "mediumpurple3")) +
  scale_y_continuous(breaks = scales::pretty_breaks(3)) +
  facet_wrap(~gene_name, scales = "free", ncol = 3) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank()) +
  labs(x = NULL, y = "TPM", fill = NULL)

ggsave("./plots/ig_c_genes.png", width = 6, height = 4.2)



# Comparison between log2 FC

condition_df <- logfc %>%
    separate(condition_id, c("time", "stim"), sep = "_") %>%
    select(-cpm) %>%
    pivot_wider(names_from = time, values_from = log2fc)

stim_diffs <- condition_df %>%
    filter(`24hr` < -60 | `72hr` < -60) %>%
    left_join(gene_names)


ggplot(anti_join(condition_df, stim_diffs),
       aes(`24hr`, `72hr`)) +
    geom_abline(linetype = 2, color = "grey") +
    geom_hline(yintercept = 0, linetype = 2, color = "grey") +
    geom_vline(xintercept = 0, linetype = 2, color = "grey") +
    geom_point(alpha = .25) +
    geom_point(data = stim_diffs,
                    aes(`24hr`, `72hr`, color = gene_type),
               alpha = .5) +
    scale_color_gdocs() +
    facet_wrap(~stim) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(color = "Gene type") +
    guides(color = guide_legend(override.aes = list(alpha = 1)))

ggsave("./plots/fc_24vs72.png", width = 6, height = 2.5)


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
    scale_color_manual(values = c("16hr_resting" = "grey70",
                                  "24hr_IgG" = "gold2",
                                  "72hr_IgG" = "gold3",
                                  "24hr_RSQ" = "mediumpurple1",
                                  "72hr_RSQ" = "mediumpurple3")) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(color = NULL)

ggsave("./plots/pca_bcell_expression.png", height = 3, width = 5)


