library(tidyverse)
library(readxl)
library(hlaseqlib)
library(ggbeeswarm)

quant <- read_tsv("./geuvadis_salmon_quants_ebv.bed")

ebv_quant <- quant %>%
    filter(chr == "ebv") %>%
    select(gid, starts_with("ERR188")) %>%
    pivot_longer(-gid, names_to = "ena_id", values_to = "tpm") %>%
    left_join(geuvadis_info, by = "ena_id") %>%
    select(sampleid = name, lab, gene_id = gid, tpm)

ebv_copies <- read_excel("./ebv_copynumbers_pone.0179446.xlsx")

ebv_df <- left_join(ebv_quant, ebv_copies, by = c("sampleid" = "samples")) %>%
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


ggplot(distinct(ebv_df, sampleid, lab, .keep_all = TRUE),
       aes(lab, `EBV load`)) +
    geom_jitter() +
    theme_bw()


# HLA

hla_gts <- polypheme_pag %>%
    mutate(allele = sub("^([^/]+).*$", "\\1", allele),
           allele = hla_trimnames(allele, 2)) %>%
    inner_join(ebv_copies, by = c("subject" = "samples"))

ggplot(hla_gts, aes(reorder(allele, `EBV load`), `EBV load`)) +
    geom_quasirandom(method = "smiley", alpha = .5) +
    facet_wrap(~locus, ncol = 1, scales = "free") +
    theme_bw()
