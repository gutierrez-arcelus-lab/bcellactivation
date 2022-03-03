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


# Fold changes
# logfc <- read_tsv("./data/logfc.tsv") %>%
#     mutate(gene_type = ifelse(gene_type %in% biotypes_order, gene_type, "other"),
#            gene_type = factor(gene_type, levels = biotypes_order))
# 
# ggplot(logfc, aes(cpm, log2fc)) +
#     geom_point(alpha = .1, size = .5) +
#     facet_wrap(~condition_id, scales = "free", ncol = 2) +
#     theme_bw() +
#     theme(panel.grid = element_blank())  +
#     labs(x = "Counts per million in condition",
#          y = expression(paste("Log"[2], FC)))
# 
# ggsave("./plots/foldchange.png", width = 5)
# 
# logfc %>%
#     filter(cpm > 10) %>%
#     ggplot(aes(cpm, log2fc)) +
#     geom_point(alpha = .1, size = .5) +
#     facet_wrap(~condition_id, scales = "free", ncol = 2) +
#     theme_bw() +
#     theme(panel.grid = element_blank()) +
#     labs(x = "Counts per million in condition",
#      y = expression(paste("Log"[2], FC)),
#      title = "CPM > 10")
# 
# ggsave("./plots/foldchange_subset.png", width = 5)
# 
# logfc %>%
#     group_by(condition_id) %>%
#     summarise(downregulated = mean(log2fc < 0),
#               upregulated = mean(log2fc > 0)) %>%
#     pivot_longer(downregulated:upregulated, names_to = "direction") %>%
#     ggplot(aes(condition_id, value, fill = direction)) +
#     geom_col(position = "dodge") +
#     scale_y_continuous(labels = scales::percent) +
#     scale_fill_npg() +
#     theme_bw() +
#     theme(panel.grid = element_blank())
#     labs(x = NULL, y = "% of genes", fill = NULL)
# 
# ggsave("./plots/foldchange_summary.png", width = 5)
# 
# # Correlations between conditions
# 
# gene_v2_df <- gene_df %>%
#     mutate(log2_cpm = log2(cpm + 1)) %>%
#     select(-tpm, -cpm) %>%
#     pivot_wider(names_from = condition_id, values_from = log2_cpm) %>%
#     pivot_longer(-(1:2), names_to = "condition", values_to = "log2_cpm") %>%
#     rename(log2_resting = `16hr_resting`) %>%
#     mutate(fc = log2_cpm - log2_resting)
#     
# cor_df <- gene_v2_df %>%
#     group_by(condition) %>%
#     summarise(y = max(log2_cpm),
#               x = min(log2_resting),
#               rho = cor(log2_cpm, log2_resting, method = "spearman"),
#               rho = round(rho, 2),
#               rho_label = paste("rho == ", rho)) %>%
#     ungroup()
# 
# fc_plot_1 <- ggplot(gene_v2_df, aes(log2_resting, log2_cpm)) +
#     geom_abline() +
#     geom_point(alpha = .1, size = .25) +
#     geom_text(data = cor_df, aes(x, y * 1.1, label = rho_label),
#              hjust = "inward", vjust = "inward",
#              parse = TRUE, size = 3) +
#     facet_wrap(~condition, scales = "free", ncol = 2) +
#     theme_bw() +
#     theme(text = element_text(size = 9),
#           plot.title = element_text(size = 9),
#           axis.title = element_text(size = 9)) +
#     labs(title = " ",
#          x = expression(paste("Log"[2], "CPM + 1 (Resting)")),
#          y = expression(paste("Log"[2], "CPM + 1")))
# 
# fc_plot_2 <- gene_v2_df %>%
#     group_by(condition) %>%
#     summarise(n = sum(abs(fc) > 2L)/n()) %>%
#     mutate(condition = factor(condition, 
#                               levels = c("24hr_IgG", "72hr_IgG", 
#                                          "24hr_RSQ", "72hr_RSQ"))) %>%
#     ggplot(aes(n, condition)) +
#     geom_col() +
#     scale_x_continuous(labels = scales::percent_format(accuracy = .1)) +
#     theme_bw() +
#     theme(text = element_text(size = 9),
#           plot.title = element_text(size = 9),
#           axis.title = element_text(size = 9)) +
#     labs(title = expression(paste("% of genes with Log"[2], " FC>2")), 
#          x = " ",
#          y = NULL)
# 
# plot_grid(fc_plot_1, NULL, fc_plot_2, nrow = 1, rel_widths = c(1, .05, .6))
# 
# 
# ggsave("./plots/scatter_resting_conditions.png", width = 6, height = 3)
# 
# 

# # Plot genes in BCR or TRL7 pathways
# 
# bcr_genes <- c("PTPN22", "CSK", "BANK1", "BLK", "LYN", "IKZF1", "IKZF3")
# tlr_genes <- c("IRAK1", "TLR7", "MYD88", "UNC93B1")
# 
# 
# pathways_df <- gene_df %>%
#     left_join(gene_names) %>%
#     select(gene_name, gene_id, condition_id, tpm) %>%
#     mutate(pathway = case_when(gene_name %in% bcr_genes ~ "BCR",
#                                gene_name %in% tlr_genes ~ "TLR7",
#                                TRUE ~ "Other")) %>%
#     filter(pathway %in% c("BCR", "TLR7"))
# 
# bcr_df <- pathways_df %>%
#     filter(condition_id %in% c("16hr_resting", "24hr_IgG", "72hr_IgG"),
#            pathway == "BCR")
# 
# tlr_df <- pathways_df %>%
#     filter(condition_id %in% c("16hr_resting", "24hr_RSQ", "72hr_RSQ"),
#            pathway == "TLR7")
# 
# bcr_plot <- ggplot(bcr_df, aes(gene_name, tpm)) +
#     geom_col(aes(fill = condition_id, alpha = condition_id), 
#              fill = "midnightblue", 
#              position = "dodge") +
#     scale_alpha_manual(values = c(.33, .66, 1)) +
#     scale_y_log10() +
#     theme_bw() +
#     theme(panel.grid = element_blank()) +
#     labs(x = NULL, y = "TPM", alpha = NULL,
#          title = "BCR pathway")
#     
# tlr_plot <- ggplot(tlr_df, aes(gene_name, tpm)) +
#     geom_col(aes(fill = condition_id, alpha = condition_id), 
#              fill = "tomato3", 
#              position = "dodge") +
#     scale_alpha_manual(values = c(.33, .66, 1)) +
#     scale_y_log10() +
#     theme_bw() +
#     theme(panel.grid = element_blank()) +
#     labs(x = NULL, y = "TPM", alpha = NULL,
#          title = "TLR pathway")
# 
# pathways_plot <- plot_grid(bcr_plot, tlr_plot, ncol = 1)
# 
# ggsave("./plots/pathways.png", pathways_plot, height = 5)

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


# selected Bentham genes for seminar

selected_genes <- c("BANK1", "BLK", "CD44", "FCGR2B", "IKZF2", "JAZF1",
                    "IRAK1", "IRF5", "SOCS1", "TNFAIP3",
                    "TASL", "SLC15A4")

bentham_selected <- gene_df %>%
  left_join(gene_names) %>%
  filter(gene_name %in% selected_genes) %>%
  mutate(condition_id = factor(condition_id, levels = conditions),
         gene_name = factor(gene_name, levels = selected_genes))

bentham_selected_plot <- ggplot(bentham_selected, aes(gene_name, tpm)) +
  geom_col(aes(fill = condition_id), position = "dodge",
           color = "black", size = .15) +
  scale_fill_manual(values = c("grey70", "gold2", "gold3",
                               "mediumpurple1", "mediumpurple3"),
                    labels = function(x) sub("_", "\n", x)) +
  scale_y_continuous(breaks = scales::pretty_breaks(3)) +
  facet_wrap(~gene_name, scales = "free", nrow =2 ) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "top") +
  labs(x = NULL, y = "TPM", fill = NULL)

ggsave("./plots/bentham_selected_genes.png", bentham_selected_plot, height = 3, width = 8)







# # HLA
# hla_df <- gene_df %>%
#     left_join(gene_names) %>%
#     filter(gene_name %in% paste0("HLA-", c("A", "B", "C", "DRA", "DRB1", "DQA1",
#                                            "DQB1", "DPA1", "DPB1"))) %>%
#     mutate(condition_id = factor(condition_id, levels = conditions))
# 
# hla_plot_1 <- ggplot(hla_df, aes(gene_name, tpm)) +
#     geom_col(aes(fill = condition_id), position = "dodge") +
#     scale_fill_manual(values = c("black", "cornflowerblue", "blue",
#                                  "salmon", "tomato3"),
#                       labels = function(x) sub("_", "\n", x)) +
#     scale_y_continuous(breaks = scales::pretty_breaks(3)) +
#     facet_wrap(~gene_name, scales = "free") +
#     theme_minimal() +
#     theme(panel.grid = element_blank(),
#           axis.text.x = element_blank(),
#           legend.position = "bottom",
#           legend.direction = "horizontal",
#           legend.background = element_blank(),
#           legend.box.background = element_rect(color = "black")) +
#     labs(x = NULL, y = "TPM", 
#          fill = NULL)
# 
# hla_fc <- hla_df %>%
#     mutate(gene_name = as.character(gene_name)) %>%
#     select(-cpm) %>%
#     mutate(tpm = tpm + 1L) %>%
#     pivot_wider(names_from = condition_id, values_from = tpm) %>%
#     pivot_longer(-(1:3), names_to = "condition_id", values_to = "tpm") %>%
#     select(gene_id, gene_name, condition_id, resting_tpm = 3, tpm) %>%
#     mutate(fc = log2(tpm) - log2(resting_tpm)) %>%
#     group_by(condition_id, gene_id) %>%
#     ungroup() %>%
#     mutate(condition_id = factor(condition_id, levels = conditions))
# 
# hla_plot_2 <- ggplot(hla_fc, aes(fc, reorder_within(gene_name, fc, condition_id))) +
#     geom_col(aes(fill = fc), show.legend = FALSE) +
#     facet_wrap(~condition_id, scales = "free_y") +
#     scale_y_reordered() +
#     scale_fill_gradient2() +
#     theme_minimal() +
#     theme(panel.grid = element_blank(),
#           plot.caption = element_text(hjust = 0),
#           axis.text.y = element_text(size = 7)) +
#     labs(x = expression(paste("Log"[2], FC)),
#          y = NULL)
# 
# ggsave("./plots/hla.png", hla_plot_1, height = 4)
# ggsave("./plots/hla_fc.png", hla_plot_2, height = 4)
# 
# 
# # Activation markers
# 
# markers <- c("CD69", "CD86")
# 
# markers_tpm <- gene_df %>%
#     left_join(gene_names) %>%
#     filter(gene_name %in% markers) %>%
#     mutate(condition_id = factor(condition_id, levels = conditions))
# 
# ggplot(markers_tpm, aes(gene_name, tpm)) +
#     geom_col(aes(fill = condition_id), position = "dodge") +
#     scale_fill_manual(values = c("black", "cornflowerblue", "blue",
#                                  "salmon", "tomato3")) +
#     theme_bw() +
#     theme(panel.grid = element_blank()) +
#     labs(x = NULL, y = "TPM", fill = NULL)
# 
# ggsave("./plots/activation_markers.png")


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


# IL21

il21 <- gene_names %>%
  filter(gene_name %in% c("IL21", "IL21R", "CD40", "CD40LG",
                          "CD27", "CXCR5", "ITGAX", "TBX21")) %>%
  inner_join(gene_df) %>%
  mutate(condition_id = factor(condition_id, levels = conditions))
  
ggplot(il21, aes(gene_name, tpm)) +
  geom_col(aes(fill = condition_id), position = "dodge",
           color = "black", size = .15) +
  scale_fill_manual(values = c("grey70", "gold2", "gold3",
                               "mediumpurple1", "mediumpurple3")) +
  scale_y_continuous(breaks = scales::pretty_breaks(3)) +
  facet_wrap(~gene_name, scales = "free", ncol = 2) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(x = NULL, y = "TPM", fill = NULL)

ggsave("./plots/il21.png", width = 6, height = 5)


# # Look at genes previously found by Jing
# gene_names_v2 <- gene_names %>%
#     mutate(ensembl = sub("\\.\\d+$", "", gene_id)) 
# 
# list_dir <- "/lab-share/IM-Gutierrez-e2/Public/B_cells_alignment/PRS/TPM_genes/"
# 
# list_files <- list.files(list_dir, pattern = "snp.txt$")
# 
# prs_list_df <- file.path(list_dir, list_files) %>%
#     setNames(sub("_snp\\.txt", "", list_files)) %>%
#     map_dfr(read_tsv, col_names = FALSE, .id = "condition_id") %>%
#     left_join(gene_names_v2, by = c("X7" = "ensembl")) %>%
#     mutate(pathway = case_when(grepl("IGG", condition_id) ~ "BCR",
#                                grepl("RSQ", condition_id) ~ "TLR7")) %>%
#     distinct(pathway, gene_id, gene_name) %>%
#     add_count(gene_name, pathway) %>%
#     mutate(gene_name = ifelse(n == 1, gene_name, paste(gene_id, gene_name, sep = "-")))
# 
# bcr_prs_df <- gene_df %>%
#     inner_join(prs_list_df) %>% 
#     filter(condition_id %in% c("16hr_resting", "24hr_IgG", "72hr_IgG"),
#            pathway == "BCR") %>%
#     distinct(condition_id, gene_id, gene_name, pathway, .keep_all = TRUE)
# 
# tlr_prs_df <- gene_df %>%
#     inner_join(prs_list_df) %>%
#     filter(condition_id %in% c("16hr_resting", "24hr_RSQ", "72hr_RSQ"),
#            pathway == "TLR7") %>%
#     distinct(condition_id, gene_id, gene_name, pathway, .keep_all = TRUE)
# 
# 
# bcr_prs_plot <- ggplot(bcr_prs_df, aes(tpm, gene_name)) +
#     geom_col(aes(fill = condition_id, alpha = fct_rev(condition_id)), 
#              fill = "midnightblue", 
#              position = "dodge") +
#     scale_alpha_manual(values = rev(c(.33, .66, 1)),
#                        guide = guide_legend(reverse = TRUE)) +
#     theme_bw() +
#     theme(panel.grid = element_blank()) +
#     labs(x = "TPM", y = NULL, alpha = NULL,
#          title = "BCR pathway")
# 
# tlr_prs_plot <- ggplot(tlr_prs_df, aes(tpm, gene_name)) +
#     geom_col(aes(fill = condition_id, alpha = fct_rev(condition_id)), 
#              fill = "tomato3", 
#              position = "dodge") +
#     scale_alpha_manual(values = rev(c(.33, .66, 1)),
#                        guide = guide_legend(reverse = TRUE)) +
#     theme_bw() +
#     theme(panel.grid = element_blank()) +
#     labs(x = "TPM", y = NULL, alpha = NULL,
#          title = "TLR7 pathway")
# 
# pathways_plot <- plot_grid(bcr_prs_plot, tlr_prs_plot, ncol = 1)
# 
# ggsave("./plots/pathways_prs.png", pathways_plot, height = 6)






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


