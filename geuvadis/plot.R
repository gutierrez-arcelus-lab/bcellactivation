library(tidyverse)
library(readxl)
library(hlaseqlib)
library(ggbeeswarm)
library(pheatmap)

quant <- read_tsv("./geuvadis_salmon_quants_ebv.bed")

ebv_quant <- quant %>%
    filter(`#chr` == "ebv") %>%
    select(gene_id = id, starts_with("ERR188")) %>%
    mutate(gene_id = sub("^([^_]+).+_(\\d+)$", "\\1.\\2", gene_id)) %>%
    pivot_longer(-gene_id, names_to = "ena_id", values_to = "tpm") %>%
    left_join(geuvadis_info, by = "ena_id") %>%
    select(sample_id = name, lab, gene_id, tpm)

ebv_copies <- read_excel("./ebv_copynumbers_pone.0179446.xlsx")

ebv_df <- left_join(ebv_quant, ebv_copies, by = c("sample_id" = "samples")) %>%
    filter(!is.na(`EBV load`), pop %in% c("GBR", "CEU", "TSI", "FIN")) %>%
    group_by(gene_id) %>%
    mutate(norm_exp = GenABEL::rntransform(tpm)) %>%
    ungroup()

top_15_expressed <- ebv_df %>%
  group_by(gene_id) %>%
  summarise(m = mean(tpm)) %>%
  ungroup() %>%
  slice_max(n = 15, order_by = m) %>%
  arrange(desc(m))

ebv_df %>%
  filter(gene_id %in% top_15_expressed$gene_id) %>%
  mutate(gene_id = factor(gene_id, levels = top_15_expressed$gene_id)) %>%
  ggplot(aes(gene_id, tpm)) +
  geom_violin()

gene_labels <- ebv_df %>%
  select(gene_id) %>%
  mutate(label = sub("(\\.\\d+){1,}$", "", gene_id)) %>%
  distinct(gene_id, label) %>%
  split(.$gene_id) %>%
  map_chr(~pull(., label))


ebv_genes_plot <- ggplot(ebv_df, aes(log2(`EBV load`), norm_exp)) +
    geom_point(size = .5, alpha = .5) +
    facet_wrap(~gene_id, scales = "free_y", 
               labeller = as_labeller(gene_labels)) +
    scale_y_continuous(breaks = scales::pretty_breaks(3)) +
    theme_bw() +
    theme(text = element_text(size = 8),
          strip.text = element_text(size = 6),
          panel.grid.minor = element_blank()) +
    labs(x = "Log2 (counts)", y = "Normalized expression")

ggsave("./plots/ebv_genes.png", ebv_genes_plot, width = 10, height = 7)





# HLA

hla_gts <- polypheme_pag %>%
    mutate(allele = sub("^([^/]+).*$", "\\1", allele),
           allele = hla_trimnames(allele, 2)) %>%
    inner_join(ebv_copies, by = c("subject" = "samples"))

ggplot(hla_gts, aes(reorder(allele, `EBV load`), `EBV load`)) +
    geom_quasirandom(method = "smiley", alpha = .5) +
    facet_wrap(~locus, ncol = 1, scales = "free") +
    theme_bw()


# PCA genotypes
pca <- read_delim("./VCF/allchr.1000G.pca", delim = " ") %>%
    mutate(SampleID = str_extract(SampleID, "PC\\d+$")) %>%
    rename(pc = SampleID) %>%
    filter(pc %in% paste0("PC", 1:10)) %>%
    pivot_longer(-pc, names_to = "sample_id") %>%
    pivot_wider(names_from = pc, values_from = value)

ids_df <- read_tsv("./geuvadis_ids.tsv")
pop_df <- read_tsv("./geuvadis_covariates.tsv") %>%
    left_join(ids_df, by = c("sampleid" = "ena_id")) %>%
    select(sample_id = name, pop)

pca %>%
    left_join(pop_df) %>%
    ggplot(aes(PC1, PC2, color = pop)) +
    geom_point()

# EBV trans-eQTL

# Correcting EBV expression for copy number
nominal <- "./results/qtltools/trans.adjust.hits.txt.gz" %>%
    read_delim(delim = " ", col_types = "ccdccdddddd") %>%
    select(SNP = VID, CHR = VCHR, BP = VPOS, P = APVAL) %>%
    mutate(CHR = ifelse(CHR == "X", 23, CHR),
           CHR = as.integer(CHR))

nominal_best <- "./results/qtltools/trans.adjust.best.txt.gz" %>%
  read_delim(delim = " ", col_types = "cddc",
             col_names = c("gene", "adj_p", "p", "snp"))

nominal_best %>% arrange(p)
    
fdr <- "./results/qtltools/trans.approx.fdr.txt" %>%
    read_delim(delim = " ", col_types = "ccdccdddddd", col_names = FALSE)
    

ggplot(nominal, aes(BP, P, color = factor(CHR))) +
    geom_hline(yintercept = fdr$X8, linetype = 2, alpha = .5) +
    geom_point(show.legend = FALSE) +
    scale_y_reverse() +
    scale_color_manual(values = rep(c("midnightblue", "tomato3"), 24)[1:23]) +
    facet_grid(~CHR, scales = "free_x", space = "free", switch = "x") +
    theme_minimal() +
    theme(panel.spacing = unit(0, "lines"),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid = element_blank()) +
    labs(y = "Adjusted p-value")

ggsave("./plots/trans_manhattan.png", height = 3)


vcf <- "./results/qtltools/signif_var.vcf" %>%
    read_tsv(comment = "##") %>%
    select(snp_id = ID, starts_with("HG"), starts_with("NA")) %>%
    pivot_longer(-snp_id, names_to = "sample_id", values_to = "gt") %>%
    separate_rows(gt, sep = "\\|") %>%
    mutate(gt = as.numeric(gt)) %>%
    group_by(snp_id, sample_id) %>%
    summarise(dosage = sum(gt)) %>%
    ungroup()

ebv_expression <- read_tsv("./ebv_phenotypes_corrected.bed.gz") %>%
    filter(id == fdr$X1) %>%
    select(gene_id = id, starts_with("HG"), starts_with("NA")) %>%
    pivot_longer(-gene_id, names_to = "sample_id", values_to = "tpm")
    
qtl_df <- inner_join(vcf, ebv_expression)

ggplot(qtl_df, aes(factor(dosage), tpm)) +
    geom_violin(fill = "grey80", alpha = .5, color = "grey80") +
    geom_quasirandom(method = "smiley") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(x = "Dosage", y = "STD normal corrected expression")
    
ggsave("./plots/trans_qtl.png", height = 4, width = 6)


# Not correcting EBV expression for copy-number
fdr_2 <- "./results/qtltools/trans.approx.fdr.v2.txt" %>%
  read_delim(delim = " ", col_types = "ccdccdddddd", col_names = FALSE)

vcf_2 <- "./results/qtltools/signif_var_v2.vcf" %>%
  read_tsv(comment = "##") %>%
  select(snp_id = ID, starts_with("HG"), starts_with("NA")) %>%
  pivot_longer(-snp_id, names_to = "sample_id", values_to = "gt") %>%
  filter(snp_id %in% fdr_2$X4) %>%
  separate_rows(gt, sep = "\\|") %>%
  mutate(gt = as.numeric(gt)) %>%
  group_by(snp_id, sample_id) %>%
  summarise(dosage = sum(gt)) %>%
  ungroup()

ebv_expression_2 <- read_tsv("./ebv_phenotypes_corrected_v2.bed.gz") %>%
  filter(id == fdr$X1) %>%
  select(gene_id = id, starts_with("HG"), starts_with("NA")) %>%
  pivot_longer(-gene_id, names_to = "sample_id", values_to = "tpm")

qtl_df_2 <- inner_join(vcf_2, ebv_expression_2)

ggplot(qtl_df_2, aes(factor(dosage), tpm)) +
  #geom_violin(fill = "grey80", alpha = .5, color = "grey80") +
  geom_quasirandom(method = "smiley", size = .5) +
  facet_wrap(~snp_id) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "Dosage", y = "STD normal corrected expression",
       title= "BALF2 gene")

ggsave("./plots/trans_qtl_2.png", height = 4, width = 6)



# Correlations
human_gene_types <- read_tsv("../data/gencode_v38_genetypes.tsv")

human_ebv_cor <- read_tsv("./human_ebv_correlations_tpm5in75.tsv") %>% 
    mutate(bonferr = p < 0.01/nrow(.))

top_50 <- human_ebv_cor %>%
    filter(abs(estimate) > .6, bonferr == TRUE) %>%
    arrange(p, desc(abs(estimate))) %>%
    select(estimate, p, gene_id_human, gene_name_human, gene_name_ebv, protein_name_ebv) %>%
    left_join(human_gene_types, 
              by = c("gene_id_human" = "gene_id", "gene_name_human" = "gene_name")) %>%
    mutate(gene_type = ifelse(gene_type == "protein_coding", gene_type, "non_protein_coding"),
           direction = ifelse(estimate > 0, "Positive", "Negative")) %>%
    mutate(lab = paste(gene_name_human, gene_name_ebv, sep = " x "),
           lab = fct_inorder(lab)) %>%
    slice_max(n = 50, order_by = estimate)


ebv_volcano <- ggplot(human_ebv_cor, aes(estimate, -log10(p))) +
    geom_point(aes(color = bonferr), size = .5, alpha = .5, show.legend = FALSE) +
    scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "black")) +
    theme_bw() +
    labs(x = "Pearson correlation coefficient")

ggsave("./plots/ebv_volcano.png", ebv_volcano, width = 4, height = 3)

ggplot(top_50, aes(estimate, lab)) +
    geom_point(aes(color = direction), size = 2, show.legend = FALSE) +
    scale_x_continuous(labels = function(x) round(x, 2)) +
    scale_color_manual(values = c("Positive" = "tomato3", "Negative" = "cornflowerblue")) +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
          axis.text.y = element_text(size = 6)) +
    labs(x = "Pearson correlation coefficient", y = NULL)

ggsave("./plots/ebv_pairs_corr.png", width = 6, height = 5)

cor_bonferr_mat <- human_ebv_cor %>%
    unite(ebv, c("gene_name_ebv", "protein_name_ebv"), sep = ": ") %>%
    unite(human, c("gene_id_human", "gene_name_human"), sep = ": ") %>%
    select(human, ebv, estimate) %>%
    pivot_wider(names_from = ebv, values_from = estimate) %>%
    as.data.frame() %>%
    column_to_rownames("human") %>%
    as.matrix()

png("./plots/ebv_heatmap.png", width = 12, height = 16, res = 300, units = "in")
pheatmap(cor_bonferr_mat, show_rownames = FALSE)
dev.off()


ebna2 <- human_ebv_cor %>%
  filter(bonferr == TRUE, gene_name_ebv == "EBNA-2") %>%
  arrange(p, desc(abs(estimate))) %>%
  select(estimate, p, gene_id_human, gene_name_human, gene_name_ebv, protein_name_ebv) %>%
  left_join(human_gene_types, 
            by = c("gene_id_human" = "gene_id", "gene_name_human" = "gene_name")) %>%
  mutate(gene_type = ifelse(gene_type == "protein_coding", gene_type, "non_protein_coding"),
         direction = ifelse(estimate > 0, "Positive", "Negative")) %>%
  slice_max(n = 50, order_by = estimate) %>%
  arrange(desc(estimate)) %>%
  mutate(gene_name_human = fct_inorder(gene_name_human))

ggplot(ebna2, aes(abs(estimate), gene_name_human)) +
  geom_point(aes(color = direction), size = 2, show.legend = FALSE) +
  scale_x_continuous(labels = function(x) round(x, 2)) +
  scale_color_manual(values = c("Positive" = "tomato3", "Negative" = "cornflowerblue")) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        axis.text.y = element_text(size = 6)) +
  labs(x = "Pearson correlation coefficient", y = NULL)

ggsave("./plots/ebna2_pairs_corr.png", width = 6, height = 4)


# Copy number

cor_copy <- read_tsv("./human_expression_ebvcopies_corr.tsv") %>%
  mutate(bonferr = p < 0.05/nrow(.))


cor_copy %>%
  filter(bonferr == TRUE) %>%
  arrange(desc(abs(estimate))) %>%
  mutate(direction = ifelse(estimate < 0, "Negative", "Positive")) %>%
  slice_max(n = 50, order_by = abs(estimate)) %>%
  arrange(direction, desc(abs(estimate))) %>%
  mutate(gid = fct_inorder(gid)) %>%
  ggplot(aes(estimate, gid)) +
  geom_col(aes(fill = direction), show.legend = FALSE) +
  scale_x_continuous(limits = c(-.4, .4)) +
  scale_fill_manual(values = c("Positive" = "tomato3", "Negative" = "cornflowerblue")) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        axis.text.y = element_text(size = 6)) +
  labs(x = "Pearson correlation coefficient", y = NULL)
  
ggsave("./plots/copies_pairs_corr.png", width = 6, height = 4)


gwas_genes <- read_lines("../bcell_scrna/reported_genes.tsv")

cor_copy %>%
  filter(gid %in% gwas_genes) %>%
  arrange(desc(abs(estimate))) %>%
  mutate(direction = ifelse(estimate < 0, "Negative", "Positive")) %>%
  slice_max(n = 50, order_by = abs(estimate)) %>%
  arrange(direction, desc(abs(estimate))) %>%
  mutate(gid = fct_inorder(gid)) %>%
  ggplot(aes(estimate, gid)) +
  geom_col(aes(fill = -log10(p))) +
  scale_x_continuous(limits = c(-.35, .35)) +
  scale_fill_viridis_c(option = "magma", limits = c(0, 10)) +
  #scale_fill_manual(values = c("Positive" = "tomato3", "Negative" = "cornflowerblue")) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        axis.text.y = element_text(size = 6)) +
  labs(x = "Pearson correlation coefficient", y = NULL)

ggsave("./plots/corr_ebv_slegenes.png", width = 6, height = 4)


## test

read_delim("./interaction/phenotypes_pcs.pca_stats", delim = " ", 
           col_names = FALSE, n_max = 3) %>%
  filter(X1 == "prop_var") %>%
  pivot_longer(-1, names_to = "pc") %>%
  mutate(pc = parse_number(pc) - 1) %>%
  ggplot(aes(pc, value)) +
  geom_line() +
  scale_x_continuous(limits = c(0, 20))



#CR2
# cr2_genos <- read_tsv("./cr2.vcf", comment = "##") %>%
#     pivot_longer(-(1:9), names_to = "sample_id", values_to = "gt") %>%
#     unite("snp_id", c("#CHROM", "POS", "REF", "ALT"), sep = "_") %>%
#     select(sample_id, snp_id, gt) %>%
#     separate_rows(gt, sep = "\\|", convert = TRUE) %>%
#     group_by(sample_id, snp_id) %>%
#     summarise(dose = sum(gt)) %>%
#     ungroup()
# 
# cr2_plot_df <- inner_join(cr2_genos, ebv_copies, by = c("sample_id" = "samples"))
#   
# cor_df <- lm(log2(`EBV load`) ~ dose, data = cr2_plot_df) %>%
#     summary() %>%
#     broom::tidy() %>%
#     filter(term == "dose") %>%
#     select(b = estimate, p.value) %>%
#     mutate(lab = paste("beta ==", round(b, 2), "*','~pvalue ==", round(p.value, 2)))
# 
# 
# ggplot(cr2_plot_df, aes(factor(dose), log2(`EBV load`))) +
#     geom_violin(fill = "grey95", color = "grey70") +
#     geom_quasirandom(method = "smiley", alpha = .5) +
#     stat_summary(fun = "median", geom = "point", shape = "\U2014", 
#                  size = 20, color = "tomato3") +
#     geom_text(data = cor_df, 
#               aes(x = 0, y = max(log2(cr2_plot_df[["EBV load"]])), label = lab), 
#               parse = TRUE, hjust = "inward", vjust = "inward", 
#               nudge_x = .5, nudge_y = .5) +
#     theme_bw() +
#     labs(x = paste("Genotypes at", paste0("chr", unique(cr2_plot_df$snp_id))),
#          y = expression(paste("log"[2], "(EBV load)")))
# 
# ggsave("./plots/cr2_sQTL.png", width = 5, height = 4)
# 
# cr2_quants <- read_tsv("./results/salmon/cr2_quants.tsv") %>%
#     left_join(select(geuvadis_info, sampleid = ena_id, name)) %>%
#     select(sample_id = name, transcript_id, tpm) %>%
#     arrange(sample_id, transcript_id)
# 
# cr2_plot_df2 <- inner_join(cr2_quants, cr2_plot_df)
# 
# cor_df2 <- cr2_plot_df2 %>%
#     mutate(tpm = log2(tpm),
#            tpm = ifelse(tpm == -Inf | tpm == Inf, 0, tpm),
#            ebv = log2(`EBV load`)) %>%
#     group_by(transcript_id) %>%
#     summarise(r = round(cor(tpm, ebv), 2)) %>%
#     ungroup()
# 
# ggplot(cr2_plot_df2, aes(log2(`EBV load`), tpm)) +
#     geom_point(size = .7) +
#     scale_y_continuous(trans = scales::pseudo_log_trans(base = 2)) +
#     facet_wrap(~transcript_id, scales = "free") +
#     theme_bw() +
#     labs(x = expression(paste("log"[2], "(EBV load)")),
#          y = expression(paste("log"[2], "(TPM)")))
# 
# ggsave("./plots/cr2_isoforms_ebv.png", width = 5, height = 4)
