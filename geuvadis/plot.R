library(tidyverse)
library(readxl)
library(hlaseqlib)
library(ggbeeswarm)

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
    filter(!is.na(`EBV load`), pop %in% c("GBR", "CEU", "TSI", "FIN"))


ebv_genes_plot <- ggplot(ebv_df, aes(`EBV load`, tpm)) +
    geom_point(size = .5, alpha = .5) +
    facet_wrap(~gene_id, scales = "free_y") +
    scale_x_continuous(breaks = scales::pretty_breaks(3)) +
    scale_y_continuous(breaks = scales::pretty_breaks(3)) +
    theme_bw() +
    theme(text = element_text(size = 8),
          panel.grid.minor = element_blank()) +
    labs(y = "TPM")

ggsave("./plots/ebv_genes.png", ebv_genes_plot, width = 10)
#ggsave("./plots/ebv_genes_log.png", ebv_genes_plot, width = 10, height = 7)



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

library(qqman)


nominal <- "./results/qtltools/trans.adjust.hits.txt.gz" %>%
    read_delim(delim = " ", col_types = "ccdccdddddd") %>%
    select(SNP = VID, CHR = VCHR, BP = VPOS, P = APVAL) %>%
    mutate(CHR = ifelse(CHR == "X", 23, CHR),
           CHR = as.integer(CHR))
    
    
fdr <- "./results/qtltools/trans.approx.fdr.txt" %>%
    read_delim(delim = " ", col_types = "ccdccdddddd", col_names = FALSE)
    
png("./plots/trans_manhattan.png", width = 8, height = 3, units = "in", res = 200)
manhattan(nominal, suggestiveline = FALSE, genomewideline = FALSE,
          col = c("blue4", "orange3"), highlight = fdr$X4)
dev.off()


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




