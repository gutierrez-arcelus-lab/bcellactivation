library(tidyverse)
library(rvest)


#plot colors

condition_colors <- c("16hr_resting" = "grey60", 
                      "24hr_IgG" = "gold2",
                      "72hr_IgG" = "gold3",
                      "24hr_RSQ" = "mediumpurple1",
                      "72hr_RSQ" = "mediumpurple3")

conditions <- names(condition_colors)

ase_df <- read_tsv("./results/ase/ase_compiled.tsv") %>%
    filter(method %in% c("ASEReadCounter", "QTLtools.75")) %>%
    mutate(method = recode(method, "QTLtools.75" = "QTLtools"),
           id = sub("^20210615_(\\d+hr_[^_]+)_.*$", "\\1", id),
           id = factor(id, levels = conditions),
           eff_size = abs(0.5 - ref_n/total))

ase_summ <- ase_df %>%
    mutate(bin = case_when(total >= 10 & total < 20 ~ 10,
                           total >= 20 & total < 30 ~ 20,
                           total >= 30 & total < 40 ~ 30,
                           total >= 40 & total < 50 ~ 40,
                           total >= 50 & total < 100 ~ 50,
                           total >= 100 & total < 1000 ~ 100,
                           total >= 1000 ~ 1000)) %>%
    group_by(method, id, bin) %>%
    summarise(n = n(),
              `No cutoff` = mean(q < 0.05),
              `AI > 0.2` = mean(q < 0.05 & eff_size > 0.2)) %>%
    ungroup() %>%
    pivot_longer(-(1:4), names_to = "cutoff", values_to = "prop")

ggplot(ase_summ, aes(bin, prop, color = id, linetype = cutoff)) +
    geom_point() +
    geom_line() +
    scale_color_manual(values = condition_colors) +
    scale_x_log10(breaks = c(10, 20, 30, 40, 50, 100, 1000)) +
    scale_y_continuous(labels = scales::percent) +
    facet_wrap(~method) +
    guides(color = FALSE) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    labs(x = "Read depth", 
         y = "% Significant (FDR < 5%)",
         linetype = "Imbalance\ncutoff")

ggsave("./plots/imbalance_cutoffs.png", width = 8, height = 3)


# Manhattan-like plot

ase_manhat_df <- ase_df %>%
    filter(method == "QTLtools") %>%
    mutate(ref_r = ref_n/total) %>%
    select(id, chr, pos, ref_r, q) %>%
    mutate(id = sub("^20210615_(\\d+hr_[^_]+)_.*$", "\\1", id),
           chr = sub("chr", "", chr),
           chr = factor(chr, levels = c(1:22, "X")))

manhat_ix <- ase_manhat_df %>%
    distinct(chr, pos) %>%
    arrange(chr, pos) %>%
    group_by(chr) %>%
    mutate(ix = 1:n()) %>%
    ungroup() 

ase_manhat_df_ix <- left_join(ase_manhat_df, manhat_ix, by = c("chr", "pos"))


ggplot(ase_manhat_df_ix, aes(ix, ref_r, color = q < 0.05)) +
    geom_point(size = .2, alpha = .5, show.legend = FALSE) +
    scale_y_continuous(breaks = c(0, .5, 1)) +
    scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "black")) +
    facet_grid(id~chr, 
               scales = "free_x", 
               space = "free_x",
               switch = "x") +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_blank(),
          panel.spacing.x = unit(0, "lines"),
          strip.text.x = element_text(size = 7),
          strip.text.y = element_text(angle = 0)) +
    labs(x = NULL, y = "REF ratio", 
         caption = "Significant ASE (FDR = 5%) in black")

ggsave("./plots/ref_ratio.png", width = 8, height = 3)


# SLE genes

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


genes_ase <- ase_df %>%
    filter(method == "QTLtools", !is.na(annot)) %>%
    separate_rows(annot, sep = ";") %>%
    separate(annot, c("gene_id", "gene_name"), sep = ":") %>%
    filter(gene_name %in% bentham_genes) %>%
    mutate(ref_r = ref_n/total) %>%
    group_by(chr, pos, ref, alt) %>%
    filter(all(conditions %in% id), any(eff_size > .2)) %>%
    ungroup() %>%
    unite("var_id", c("gene_name", "pos", "ref", "alt")) %>%
    select(id, var_id, ref_r)


ggplot(genes_ase, aes(id, ref_r, fill = id)) +
    geom_hline(yintercept = .5, linetype = 2) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = condition_colors) +
    scale_y_continuous(breaks = c(0, .5, 1)) +
    facet_wrap(~var_id, scales = "free_x", ncol = 5) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.text = element_text(size = 6)) +
    labs(x = NULL, y = "REF allele ratio", fill = NULL)

ggsave("./plots/sle_genes_ai.png", height = 8, width = 8)













