library(clusterProfiler)
library(org.Hs.eg.db)
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



# SLE genes
gwas_genes <- read_tsv("../bcell_scrna/reported_genes.tsv") %>%
    mutate(gene = recode(gene, "CXorf21" = "TASL"))

genes_ase <- ase_df %>%
    filter(method == "ASEReadCounter", !is.na(annot)) %>%
    select(-method) %>%
    separate_rows(annot, sep = ";") %>%
    separate(annot, c("gene_id", "gene_name"), sep = ":") %>%
    mutate(chr = sub("chr", "", chr),
           var_id = paste0(gene_name, ":", chr, "_", pos, "_", ref, "_", alt)) %>%
    distinct(stim = id, var_id, gene_name, gene_id, pos, ref_ratio, q)

ase_igg <- genes_ase %>%
    filter(stim %in% c("16hr_resting", "24hr_IgG", "72hr_IgG")) %>%
    arrange(var_id) %>%
    group_by(var_id) %>%
    filter(all(c("16hr_resting", "24hr_IgG", "72hr_IgG") %in% stim)) %>%
    mutate(eff = abs(first(ref_ratio) - ref_ratio)) %>%
    filter(any(eff > .2 & q < .05)) %>%
    ungroup() %>%
    select(-eff)

ase_rsq <- genes_ase %>%
    filter(stim %in% c("16hr_resting", "24hr_RSQ", "72hr_RSQ")) %>%
    arrange(var_id) %>%
    group_by(var_id) %>%
    filter(all(c("16hr_resting", "24hr_RSQ", "72hr_RSQ") %in% stim)) %>%
    mutate(eff = abs(first(ref_ratio) - ref_ratio)) %>%
    filter(any(eff > .2 & q < .05)) %>%
    ungroup() %>%
    select(-eff)

bind_rows("BCR stim" = ase_igg, "TLR7 stim" = ase_rsq, .id = "condition") %>%
    filter(gene_name %in% gwas_genes$gene) %>%
    ggplot(aes(stim, var_id)) +
    geom_tile(aes(fill = ref_ratio)) +
    scale_fill_gradient2(midpoint = 0.5) +
    facet_wrap(~condition, scales = "free") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          strip.text = element_text(size = 12),
          strip.background = element_blank()) +
    labs(fill = "Ref allele\nratio")

ggsave("./plots/ase_gwasgenes.png", width = 8, height = 5, dpi = 300)

sle_eur_genes <- gwas_genes %>%
    filter(grepl("Langefeld|Bentham", author)) %>%
    distinct(gene_name = gene) %>%
    inner_join(genes_ase) %>%
    filter(q < 0.05) %>%
    count(gene_name, stim) 

sle_eur_genes %>% print(n = Inf)

sle_eur_genes %>%
    count(stim, sort = TRUE)


sle_eur_genes %>%
    complete(gene_name, stim, fill = list(n = 0)) %>%
    mutate(stim = recode(stim, "16hr_resting" = "Resting", 
                         "24hr_IgG" = "BCR 24hr",
                         "72hr_IgG" = "BCR 72hr",
                         "24hr_RSQ" = "TLR7 24hr",
                         "72hr_RSQ" = "TLR7 72hr")) %>%
    arrange(gene_name, stim) %>%
    ggplot(aes(stim, gene_name)) +
    geom_point(aes(color = n), size = 4, alpha = .25) +
    geom_text(aes(label = n), size = 3) +
    scale_color_viridis_c(guide = guide_colorbar(barwidth = .5)) +
    theme_minimal() +
    theme(axis.title = element_blank(),
          panel.grid = element_line(color = "grey96"),
          plot.background = element_rect(fill = "white"),
          plot.title = element_text(size = 8)) +
    labs(title = "Number of ASE variants per condition for SLE genes")

ggsave("ase_lupus.png", width = 5)



## enrichment
run_enrichment <- function(gene_list) {
    
    enrichGO(gene = gene_list,
             OrgDb = org.Hs.eg.db,
             ont = "BP",
             keyType = "ENSEMBL",
             pAdjustMethod = "fdr",
             pvalueCutoff = 0.05,
             qvalueCutoff = 0.05,
             readable = TRUE) %>%
        as.data.frame() %>%
        as_tibble()
}

go_igg <- ase_igg %>%
    pull(gene_id) %>%
    unique() %>%
    str_remove("\\.\\d+$") %>%
    run_enrichment()

go_rsq <- ase_rsq %>%
    pull(gene_id) %>%
    unique() %>%
    str_remove("\\.\\d+$") %>%
    run_enrichment()

bind_rows("BCR" = go_igg, "TLR7" = go_rsq, .id = "condition") %>%
    select(condition, Description, Count, qvalue) %>%
    ggplot(aes(condition, reorder(Description, -log10(qvalue)))) +
    geom_point(aes(size = Count, fill = -log10(qvalue)), shape = 21) +
    scale_fill_viridis_c(option = "magma") +
    theme_bw() +
    theme(axis.title = element_blank()) +
    labs(size = "# of Genes:")

ggsave("./plots/ase_go.png", dpi = 300)



##
genes_ase %>%
    mutate(gene_name = sub("^([^_]+).*$", "\\1", var_id)) %>%
    select(condition = id, gene_name, var_id, qvalue = q) %>%
    mutate(condition = factor(condition, levels = conditions)) %>%
    arrange(gene_name, var_id, condition) %>%
    write_tsv("./results/ase/ase_pvalues_SLEgenes.tsv")
    
## 
all_genes <- ase_df %>%
    filter(method == "ASEReadCounter", !is.na(annot)) %>%
    separate_rows(annot, sep = ";") %>%
    separate(annot, c("gene_id", "gene_name"), sep = ":") %>%
    select(-method) %>%
    arrange(gene_name, chr, pos, id) %>%
    group_by(chr, pos, ref, alt) %>%
    filter(all(conditions %in% id) &
               first(eff_size < .1) &
               any(depth > 20 & eff_size > .2 & q < 0.05)) %>%
    ungroup() %>%
    unite("var_id", c("gene_name", "pos", "ref", "alt")) 

all_genes %>%
    mutate(score = eff_size - first(eff_size)) %>%
    group_by(var_id) %>%
    slice(which.max(score)) %>%
    ungroup() %>%
    select(var_id, score) %>%
    arrange(desc(score)) %>%
    print(n = 30)

selected_vars <- ase_df %>%
    filter(method == "ASEReadCounter") %>%
    filter((chr == "chr20" & pos == 25447658) |
           (chr == "chr8" & pos == 11494547) |
           (chr == "chr19" & pos == 7696361) |
           (chr == "chr19" & pos == 10351440)) %>%
    separate_rows(annot, sep = ";") %>%
    separate(annot, c("gene_id", "gene_name"), sep = ":") %>%
    select(-method) %>%
    arrange(gene_name, chr, pos, id) %>%
    mutate(gene_name = ifelse(is.na(gene_name), "FCER2", gene_name)) %>%
    unite("var_id", c("gene_name", "pos", "ref", "alt")) %>%
    mutate(alt_n = depth - ref_n) %>%
    select(id, var_id, REF = ref_n, ALT = alt_n) %>%
    pivot_longer(REF:ALT, names_to = "allele", values_to = "count") %>%
    mutate(allele = factor(allele, levels = c("REF", "ALT")))

ggplot(selected_vars, aes(id, count, fill = id)) +
    geom_bar(aes(linetype = allele), 
             stat = "identity", position = "dodge", size = .5, 
             color = "black") +
    scale_fill_manual(values = condition_colors) +
    scale_y_continuous(breaks = scales::pretty_breaks(3)) +
    scale_linetype_manual(values = c("REF" = 1, "ALT" = 3)) +
    facet_wrap(~var_id, scales = "free_y") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_blank()) +
    labs(x = NULL, y = "Allele counts", fill = "stim") +
    guides(linetype = guide_legend(override.aes = list(fill = NA)),
           fill = guide_legend(override.aes = list(linetype = 0)))

ggsave("./plots/selected_vars_ai.png", width = 6, height = 3)



# ASE gene-level
ase_gene <- read_tsv("./results/phaser/MG8989_gene_ae.txt") %>%
    mutate(bam = sub("20210615_(\\d+hr_[^_]+)_.+$", "\\1", bam),
           bam = factor(bam, levels = conditions)) %>%
    extract(name, c("gene_id", "gene_name"), "(ENSG[0-9.]+)_(.+)") %>%
    select(condition = bam, everything()) %>%
    arrange(contig, start, gene_name, condition)
    
    
ase_gene_selected <- ase_gene %>%
    filter(gene_name %in% c("C22orf34", "FCER2")) %>%
    select(condition, gene_name, aCount, bCount, variants) %>%
    rename(Hap1 = aCount, Hap2 = bCount) %>%
    pivot_longer(Hap1:Hap2, names_to = "hap", values_to = "count") %>%
    mutate(variants = gsub("chr\\d+_", "", variants),
           variants = paste(gene_name, variants, sep = ": "))
    
# check same variants in the SNP-based data
ase_vars_df <- ase_gene_selected %>%
    distinct(gene_name, variants) %>%
    mutate(variants = sub("^([A-Za-z0-9]+: )", "", variants)) %>%
    separate_rows(variants, sep = ",")  %>%
    separate(variants, c("pos", "ref", "alt"), sep = "_", convert = TRUE) %>%
    inner_join(distinct(ase_gene, chr = contig, gene_name))


separate_variants <- ase_df %>% 
    filter(method == "ASEReadCounter") %>%
    select(-method) %>%
    inner_join(ase_vars_df, by = c("chr", "pos", "ref", "alt")) %>%
    unite(var_id, c("pos", "ref", "alt"), sep = "_") %>%
    mutate(alt_n = depth - ref_n) %>%
    select(condition = id, var_id, REF = ref_n, ALT = alt_n, gene_name) %>%
    pivot_longer(REF:ALT, names_to = "allele", values_to = "count") %>%
    complete(condition, 
             nesting(var_id, gene_name, allele), 
             fill = list(count = 0)) %>%
    mutate(allele = factor(allele, levels = c("REF", "ALT"))) %>%
    arrange(gene_name, condition, var_id, allele)
    





# compare them

corf_ase_plot <- ase_gene_selected %>%
    filter(gene_name == "C22orf34") %>%
    unite(label_x, c("hap", "condition"), sep = ":", remove = FALSE) %>%
    mutate(label_x = fct_inorder(label_x)) %>%
    ggplot(aes(label_x, count, fill = condition)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = condition_colors) +
    scale_x_discrete(labels = function(x) sub("^(Hap[12]):.+$", "\\1", x)) +
    scale_y_continuous(breaks = scales::pretty_breaks(3)) +
    facet_wrap(~variants, scales = "free_y", ncol = 1) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(x = NULL, y = "Haplotypic count") 

fcer2_ase_plot <- ase_gene_selected %>%
    filter(gene_name == "FCER2") %>%
    unite(label_x, c("hap", "condition"), sep = ":", remove = FALSE) %>%
    mutate(label_x = fct_inorder(label_x)) %>%
    ggplot(aes(label_x, count, fill = condition)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = condition_colors) +
    scale_x_discrete(labels = function(x) sub("^(Hap[12]):.+$", "\\1", x)) +
    scale_y_continuous(breaks = scales::pretty_breaks(3)) +
    facet_wrap(~variants, scales = "free_y", ncol = 1) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(x = NULL, y = "Haplotypic count") 


corf_snp_plot <- separate_variants %>%
    filter(gene_name == "C22orf34") %>%
    unite(label_x, c("allele", "condition"), sep = ":", remove = FALSE) %>%
    mutate(label_x = fct_inorder(label_x)) %>%
    ggplot(aes(label_x, count, fill = condition)) +
    geom_bar(stat = "identity", position = "dodge", show.legend = FALSE) +
    scale_fill_manual(values = condition_colors) +
    scale_x_discrete(labels = function(x) sub("^([^:]+):.+$", "\\1", x)) +
    scale_y_continuous(breaks = scales::pretty_breaks(3)) +
    facet_wrap(~var_id) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    labs(x = NULL, y = "Allele counts")

fcer2_snp_plot <- separate_variants %>%
    filter(gene_name == "FCER2") %>%
    unite(label_x, c("allele", "condition"), sep = ":", remove = FALSE) %>%
    mutate(label_x = fct_inorder(label_x)) %>%
    ggplot(aes(label_x, count, fill = condition)) +
    geom_bar(stat = "identity", position = "dodge", show.legend = FALSE) +
    scale_fill_manual(values = condition_colors) +
    scale_x_discrete(labels = function(x) sub("^([^:]+):.+$", "\\1", x)) +
    scale_y_continuous(breaks = scales::pretty_breaks(3)) +
    facet_wrap(~var_id, ncol = 4) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(size = 7, angle = 45, hjust = 1, vjust = 1)) +
    labs(x = NULL, y = "Allele counts")


plot_grid(corf_ase_plot + theme(legend.position = "none") + labs(title = "Haplotype"), 
          NULL,
          corf_snp_plot + labs(title = "SNP-level"), 
          rel_heights = c(1, .1, 1),
          ncol = 1) %>%
    plot_grid(get_legend(corf_ase_plot), rel_widths = c(1, .2))

ggsave("./plots/gene_level_ai_c22orf.png", height = 5, width = 7)


plot_grid(fcer2_ase_plot + 
              theme(strip.text = element_text(size = 7),
                    legend.position = "none") + 
              labs(title = "Haplotype"), 
          NULL,
          fcer2_snp_plot + labs(title = "SNP-level"), 
          rel_heights = c(.75, .1, 1),
          ncol = 1) %>%
    plot_grid(get_legend(corf_ase_plot), rel_widths = c(1, .2))
    
ggsave("./plots/gene_level_ai_fcer2.png", height = 5, width = 7)


### test
ase_df %>%
    filter(method == "ASEReadCounter", !is.na(annot)) %>%
    separate_rows(annot, sep = ";") %>%
    separate(annot, c("gene_id", "gene_name"), sep = ":") %>%
    select(-method) %>%
    arrange(gene_name, chr, pos, id) %>%
    filter(gene_name == "FCER2")

gene_exp_df <- read_tsv("../data/gene_quants.tsv")

fcer_id <- ase_gene %>%
    filter(gene_name == "FCER2") %>%
    pull(gene_id) %>%
    unique()

gene_exp_df %>%
    filter(gene_id == fcer_id)

ase_gene %>%
    select(-gw_phased) %>%
    mutate(ai = abs(0.5 - aCount/totalCount)) %>%
    group_by(gene_name) %>%
    filter(all(!is.na(variants))) %>%
    filter(any(totalCount >= 20)) %>%
    group_by(gene_name, variants) %>%
    filter(first(ai) < .1 & any(ai > .2 & totalCount >= 16)) %>%
    mutate(score = ai - first(ai)) %>%
    summarise(score = max(score)) %>%
    separate_rows(variants, sep = ",") %>%
    filter(n() >= 4) %>%
    summarise(variants = paste(variants, collapse = ","),
              score = unique(score)) %>%
    ungroup() %>%
    arrange(desc(score)) %>%
    print(n = 50)

ase_gene %>% filter(gene_name %in% "TYK2") %>%
    select(-gw_phased) %>%
    mutate(ai = round(abs(0.5 - aCount/totalCount), 2))
    
    

