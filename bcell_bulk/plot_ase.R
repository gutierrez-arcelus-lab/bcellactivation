library(tidyverse)
library(rvest)
library(ggridges)
library(cowplot)


#plot colors

condition_colors <- c("16hr_resting" = "grey60", 
                      "24hr_IgG" = "gold2",
                      "72hr_IgG" = "gold3",
                      "24hr_RSQ" = "mediumpurple1",
                      "72hr_RSQ" = "mediumpurple3")

conditions <- names(condition_colors)

ase_df <- read_tsv("./results/ase/ase_compiled.tsv") %>%
    mutate(id = factor(id, levels = conditions),
           ref_ratio = ref_n/depth,
           eff_size = abs(0.5 - ref_ratio))

ase_summ <- ase_df %>%
    mutate(bin = case_when(depth >= 10 & depth < 20 ~ 10,
                           depth >= 20 & depth < 30 ~ 20,
                           depth >= 30 & depth < 40 ~ 30,
                           depth >= 40 & depth < 50 ~ 40,
                           depth >= 50 & depth < 100 ~ 50,
                           depth >= 100 & depth < 1000 ~ 100,
                           depth >= 1000 ~ 1000)) %>%
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

# Distribution of REF ratios

ase_hist_df <- ase_df %>%
    filter(method == "ASEReadCounter") %>%
    select(id, chr, pos, ref_ratio)

ase_hist_avg <- ase_hist_df %>%
    group_by(id) %>%
    summarise(m = mean(ref_ratio))

hist_1 <- ggplot(ase_hist_df, aes(ref_ratio, id, fill = id)) +
    geom_density_ridges(stat = "binline", scale = .95, bins = 40,
                        draw_baseline = FALSE,
                        center = 0.5, color = NA) +
    scale_fill_manual(values = condition_colors) +
    theme_bw() +
    labs(x = "Reference fraction", y = NULL, fill = NULL)

hist_2 <- ggplot(ase_hist_avg, aes(id, m, fill = id)) +
    geom_col() +
    scale_fill_manual(values = condition_colors) +
    scale_x_discrete(labels = rep(" ", 5)) +
    scale_y_continuous(breaks = c(0, .25, .5)) +
    theme_bw() +
    labs(x = " ", y = "Average REF ratio", fill = NULL)
    
plot_grid(hist_1 + theme(legend.position = "none"),
          hist_2 + theme(legend.position = "none"),
          get_legend(hist_1),
          rel_widths = c(1, .33, .33), nrow = 1)

ggsave("./plots/ref_fraction_hist.png", width = 6, height = 3)


ase_df %>%
    filter(method == "ASEReadCounter") %>%
    select(id, chr, pos, depth) %>%
    ggplot(aes(log10(depth), id, fill = id)) +
    geom_density_ridges(stat = "binline", scale = .95,
                        draw_baseline = FALSE,
                        color = NA,
                        bins = 40) +
    scale_fill_manual(values = condition_colors) +
    theme_bw() +
    labs(x = "Log10 read counts", y = NULL, fill = NULL)

ggsave("./plots/coverage_hist.png", width = 5, height = 4)


ase_df %>%
    count(method, id) %>%
    ggplot(aes(id, n, fill = method)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_x_discrete(labels = function(x) sub("_", "\n", x)) +
    scale_y_continuous(labels = scales::comma) +
    scale_fill_manual(values = c("grey30", "grey")) +
    theme_bw() +
    theme(text = element_text(size = 8),
          panel.grid = element_blank()) +
    labs(x = NULL,
         y = "Number of sites with >= 10 reads", fill = NULL)

ggsave("./plots/number_of_sites.png", width = 5)



# Manhattan-like plot
# 
# ase_manhat_df <- ase_df %>%
#     filter(method == "ASEReadCounter") %>%
#     select(id, chr, pos, ref_ratio, q) %>%
#     mutate(chr = sub("chr", "", chr),
#            chr = factor(chr, levels = c(1:22, "X")))
# 
# manhat_ix <- ase_manhat_df %>%
#     distinct(chr, pos) %>%
#     arrange(chr, pos) %>%
#     group_by(chr) %>%
#     mutate(ix = 1:n()) %>%
#     ungroup() 
# 
# ase_manhat_df_ix <- left_join(ase_manhat_df, manhat_ix, by = c("chr", "pos"))
# 
# 
# ggplot(ase_manhat_df_ix, aes(ix, ref_r, color = q < 0.05)) +
#     geom_point(size = .2, alpha = .5, show.legend = FALSE) +
#     scale_y_continuous(breaks = c(0, .5, 1)) +
#     scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "black")) +
#     facet_grid(id~chr, 
#                scales = "free_x", 
#                space = "free_x",
#                switch = "x") +
#     theme_minimal() +
#     theme(panel.grid = element_blank(),
#           axis.text.x = element_blank(),
#           panel.spacing.x = unit(0, "lines"),
#           strip.text.x = element_text(size = 7),
#           strip.text.y = element_text(angle = 0)) +
#     labs(x = NULL, y = "REF ratio", 
#          caption = "Significant ASE (FDR = 5%) in black")
# 
# ggsave("./plots/ref_ratio.png", width = 8, height = 3)


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
    filter(method == "ASEReadCounter", !is.na(annot)) %>%
    separate_rows(annot, sep = ";") %>%
    separate(annot, c("gene_id", "gene_name"), sep = ":") %>%
    filter(gene_name %in% bentham_genes) %>%
    group_by(chr, pos, ref, alt) %>%
    filter(all(conditions %in% id) & any(q < 0.05)) %>%
    ungroup() %>%
    unite("var_id", c("gene_name", "pos", "ref", "alt")) %>%
    mutate(alt_n = depth - ref_n) %>%
    select(id, var_id, REF = ref_n, ALT = alt_n) %>%
    pivot_longer(REF:ALT, names_to = "allele", values_to = "count") %>%
    mutate(allele = factor(allele, levels = c("REF", "ALT")))


ggplot(genes_ase, aes(id, count, fill = allele)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("REF" = "grey30", "ALT" = "grey")) +
    scale_x_discrete(labels = function(x) sub("_", "\n", x)) +
    scale_y_continuous(breaks = scales::pretty_breaks(3)) +
    facet_wrap(~var_id, scales = "free_y", ncol = 4) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = "bottom") +
    labs(x = NULL, y = "Allele counts")

ggsave("./plots/sle_genes_refalt.png", height = 8, width = 8)







    
    
    
    
    
    
    
    
    





