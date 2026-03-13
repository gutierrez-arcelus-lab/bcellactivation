library(tidyverse)
library(glue)
library(ggrepel)
library(tximport)

# leafcutter
read_leaf <- function(stim) {
    
    sig <- 
        glue("./results/{stim}_cluster_significance.txt") |>
        read_tsv() |>
	filter(status == "Success") |>
        separate(cluster, c("chr", "cluster"), sep = ":") |>
        select(cluster, p, padj = p.adjust, genes)

    eff <- 
        glue("./results/{stim}_effect_sizes.txt") |>
        read_tsv() |>
        mutate(cluster = sub("^.*(clu.*)$", "\\1", intron)) |>
        select(cluster, logef, deltapsi)

    inner_join(sig, eff, join_by(cluster))
}

leaf_df <- 
    c("TLR7", "BCR", "DN2") |> 
    {function(x) setNames(x, x)}() |>
    map_dfr(~read_leaf(.x) |>
	    group_by(cluster) |>
	    slice_min(p) |>
	    slice_max(abs(deltapsi)) |>
	    ungroup(),
	    .id = "stim") |>
    mutate(stim = fct_inorder(stim))

plot_labels <- 
    leaf_df |>
    filter(padj <= 0.05) |>
    group_by(stim) |>
    top_n(20, -log10(p)) |>
    ungroup() |>
    mutate(genes = str_trunc(genes, 10))

leaf_plot <- 
    ggplot(leaf_df, aes(x = deltapsi, y = -log10(p))) +
    geom_point(size = .5) +
    geom_text_repel(data = plot_labels, 
		    aes(x = deltapsi, y = -log10(p), label = genes),
		    size = 2,
		    min.segment.length = 0,
		    segment.size = .25,
		    max.overlaps = Inf) +
    facet_wrap(~stim, ncol = 1) +
    theme_bw() +
    theme(panel.grid = element_blank())

ggsave("./plots/leafcutter_ds.png", leaf_plot, height = 8, width = 5)

# DS genes
leaf_df |>
    filter(!is.na(genes), padj <= 0.05) |>
    distinct(stim, genes) |>
    count(stim)


# GSEA
library(fgsea)

pathwaysgo <- gmtPathways("/lab-share/IM-Gutierrez-e2/Public/References/msigdb/c5.all.v2022.1.Hs.symbols.gmt.txt")
gobp <- keep(pathwaysgo, grepl("^GOBP", names(pathwaysgo)))

kegg <- gmtPathways("/lab-share/IM-Gutierrez-e2/Public/References/msigdb/c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt")

gsea_list <- 
    leaf_df |>
    group_by(stim, genes) |>
    slice_min(p) |>
    ungroup() |>
    distinct(stim, genes, p) |>
    mutate(p = -log10(p)) |>
    separate_rows(genes, sep = ",") |>
    filter(genes %in% unique(unlist(gobp))) |>
    select(stim, SYMBOL = genes, stat = p) |>
    arrange(stim, desc(stat)) |>
    {function(x) split(x, x$stim)}() |>
    map(~select(., -stim)) |>
    map(deframe)

gsea_res <- 
    map(gsea_list, ~fgsea(pathways = gobp, stats = ., scoreType = "pos"))

bind_rows(gsea_res, .id = "stim") |>
    as_tibble() |>
    group_by(stim) |>
    slice_max(n = 10, -log10(pval)) |>
    ungroup() |>
    arrange(stim, desc(-log10(pval))) |>
    mutate(pathway = str_remove(pathway, "^GOBP_")) |> print(n = Inf)



# Salmon
tx_to_gene <- 
    file.path("/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens",
	      "gencode.v41.primary_assembly.annotation.gtf.gz") |>
    vroom::vroom(comment = "#", col_names = FALSE) |>
    filter(X3 == "transcript") |>
    mutate(tx_id = str_extract(X9, "(?<=transcript_id\\s\")[^\"]+"),
	   gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
	   gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+")) |>
    select(tx_id, gene_id, gene_name)

meta <- 
    read_tsv("./metadata_qced.tsv") |>
    separate(sample_id, c("stim", "donor", "replic"), sep = "_", remove = FALSE) |>
    select(sample_id, stim, donor)

salmon_files <- 
    file.path("./salmon_quant", meta$sample_id, "quant.sf") |>
    setNames(meta$sample_id)



# PCA
txi <- 
    tximport(salmon_files, 
	     type = "salmon", 
	     tx2gene = select(tx_to_gene, tx_id, gene_id))

sample_table <- 
    meta |>
    column_to_rownames("sample_id")
   
dds <- 
    DESeq2::DESeqDataSetFromTximport(txi, sample_table, ~stim + donor) |>
    DESeq2::estimateSizeFactors()

pca_data <-
    DESeq2::plotPCA(DESeq2::vst(dds), intgroup = "stim", ntop = 2000, returnData = TRUE) |>
    as_tibble() |>
    select(sample_id = name, stim, PC1, PC2)


salmon_df <- 
    map_dfr(salmon_files, vroom::vroom, .id = "sample_id") |>
    left_join(tx_to_gene, join_by(Name == tx_id)) |>
    left_join(meta, join_by(sample_id)) |>
    select(sample_id, donor, stim, gene_id, gene_name, tx_id = Name, tpm = TPM)

tnf_df <- 
    salmon_df |>
    filter(gene_name == "TNFSF4") |>
    mutate_at(vars(gene_id, tx_id), ~str_remove(., "\\.\\d+$")) |>
    mutate(stim = factor(stim, levels = c("unstday0", "TLR7", "BCR", "DN2")))

tnf_plot <- 
    ggplot(tnf_df, aes(x = tpm, y = tx_id, color = stim)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(position = position_dodge(width = 0.75)) +
    scale_color_manual(values = c("gray50", "forestgreen", "blue", "tomato3")) +
    facet_wrap(~stim, nrow = 1) +
    labs(x = "Transcripts per million", y = NULL)

ggsave("./plots/tnf_expression.png", tnf_plot, width = 8, height = 2)
