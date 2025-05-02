library(tidyverse)
library(glue)
library(ggrepel)
library(tidytext)
library(tximport)
library(cowplot)
library(fgsea)

fig_colors <- 
    read_tsv("../figure_colors_v2.txt", col_names = c("stim", "timep", "color")) |>
    unite("condition", c(stim, timep), sep = " ") |>
    mutate(condition = paste0(condition, "h")) |>
    deframe()

rna_colors <- fig_colors[c("Unstim 0h", "TLR7c 24h", "BCRc 24h", "DN2c 24h")]
rna_colors[c("TLR7c 24h", "BCRc 24h", "DN2c 24h")] <- fig_colors[c("TLR7c 24h", "BCRc 48h", "DN2c 48h")] 


# Fig A
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
    read_tsv("../bcell_rnaseq/4-splicing/metadata_qced.tsv") |>
    separate(sample_id, c("stim", "donor", "replic"), sep = "_", remove = FALSE) |>
    select(sample_id, stim, donor) |>
    mutate(stim = recode(stim, 
			 "unstday0" = "Unstim 0h",
			 "TLR7" = "TLR7c 24h",
			 "BCR" = "BCRc 24h",
			 "DN2" = "DN2c 24h"),
	   stim = factor(stim, levels = names(rna_colors)))

salmon_files <- 
    file.path("../bcell_rnaseq/4-splicing/salmon_quant", meta$sample_id, "quant.sf") |>
    setNames(meta$sample_id)

# PCA
txi <- 
    tximport(salmon_files, 
	     type = "salmon", 
	     tx2gene = select(tx_to_gene, tx_id, gene_id))

sample_table <- 
    column_to_rownames(meta, "sample_id")
   
dds <- 
    DESeq2::DESeqDataSetFromTximport(txi, sample_table, ~1) |>
    DESeq2::estimateSizeFactors()

pca_data <-
    DESeq2::plotPCA(DESeq2::vst(dds), intgroup = "stim", ntop = 4000, returnData = TRUE) |>
    as_tibble() |>
    select(sample_id = name, stim, PC1, PC2)

pca_plot <- 
    ggplot(pca_data, aes(x = PC1, y = PC2)) +
    geom_hline(yintercept = 0, color = "grey80", linewidth = .25) +
    geom_vline(xintercept = 0, color = "grey80", linewidth = .25) +
    geom_point(aes(fill = stim), size = 3, shape = 21, stroke = .2) +
    scale_fill_manual(values = rna_colors) +
    theme_minimal() +
    theme(
	  axis.text = element_blank(),
	  axis.title = element_text(size = 9),
	  legend.text = element_text(size = 8),
	  legend.title = element_text(size = 9),
	  panel.grid = element_blank()) +
    guides(fill = guide_legend(override.aes = list(size = 4))) +
    labs(fill = "Stim:")

fig_a_title <- 
    ggdraw() + 
    draw_label(
	       "Standard-input RNA-seq on 16 donors:\nPCA on gene expression data separates conditions",
	       x = 0,
	       size = 9,
	       hjust = 0
	       ) +
    theme(text = element_text(size = 9),
	  plot.margin = margin(l = 1.25, unit = "lines"))

fig_a_grid <- 
    plot_grid(fig_a_title, NULL, pca_plot, 
	      ncol = 1, rel_heights = c(.15, .05, 1))




# Fig B 
read_leaf <- function(stim) {
    
    sig <- 
        glue("../bcell_rnaseq/4-splicing/results/{stim}_cluster_significance.txt") |>
        read_tsv() |>
	filter(status == "Success") |>
        separate(cluster, c("chr", "cluster"), sep = ":") |>
        select(cluster, p, padj = p.adjust, genes)

    eff <- 
        glue("../bcell_rnaseq/4-splicing/results/{stim}_effect_sizes.txt") |>
        read_tsv() |>
        mutate(cluster = sub("^.*(clu.*)$", "\\1", intron)) |>
        select(cluster, logef, deltapsi)

    inner_join(sig, eff, join_by(cluster))
}

leaf_df <- 
    c("TLR7", "BCR", "DN2") |> 
    {function(x) setNames(x, x)}() |>
    map_dfr(read_leaf, .id = "stim") |>
    mutate(stim = recode(stim, "TLR7" = "TLR7c 24h", "BCR" = "BCRc 24h", "DN2" = "DN2c 24h"),
	   stim = factor(stim, levels = names(rna_colors)[-1]))

leaf_summary <- 
    leaf_df |>
    filter(padj <= 0.05, !is.na(genes)) |>
    distinct(stim, genes) |>
    count(stim)
    
leaf_summ_plot <- 
    ggplot(leaf_summary, aes(x = stim, y = n)) +
    geom_col(aes(fill = stim), color = "black", linewidth = .25) +
    scale_y_continuous(breaks = seq(0, 6000, by = 1000)) +
    scale_fill_manual(values = rna_colors) +
    theme_minimal() +
    theme(axis.text.x = element_blank()) +
    guides(fill = "none") +
    labs(x = NULL, y = "Number of genes")

fig_b_title <- 
    ggdraw() + 
    draw_label(
	       "Number of differentially spliced\ngenes (vs. Unstimulated cells)",
	       x = 0,
	       size = 9,
	       hjust = 0
	       ) +
    theme(text = element_text(size = 9),
	  plot.margin = margin(l = 1.25, unit = "lines"))

fig_b_grid <- 
    plot_grid(fig_b_title, NULL, leaf_summ_plot, 
	      ncol = 1, rel_heights = c(.15, .05, 1))



# Fig C
volcano_data <- 
    leaf_df |>
    mutate(absd = abs(deltapsi)) |>
    group_by(stim, genes) |>
    slice_min(p) |>
    slice_max(absd) |>
    ungroup() |>
    distinct(stim, genes, p, absd)

volcano_labels <- 
    volcano_data |>
    group_by(stim) |>
    top_n(15, -log10(p)) |>
    ungroup() |>
    mutate(genes = str_trunc(genes, 10))

volcano_plot <- 
    ggplot(volcano_data, aes(x = absd, y = -log10(p))) +
    geom_point(aes(color = stim), size = .5) +
    geom_text_repel(data = volcano_labels, 
		    aes(x = absd, y = -log10(p), label = genes),
		    size = 2,
		    min.segment.length = 0,
		    segment.size = .2,
		    max.overlaps = Inf,
		    fontface = "bold") +
    scale_x_continuous(breaks = c(0, .5, 1),
		       labels = c("0", "0.5", "1")) +
    scale_color_manual(values = rna_colors) +
    facet_wrap(~stim, nrow = 1) +
    theme_minimal() +
    theme(
	  axis.text = element_text(size = 8),
	  axis.title = element_text(size = 9),
	  strip.text = element_text(size = 9),
	  panel.grid.minor.y = element_blank(),
	  panel.grid.major.y = element_blank()) +
    guides(color = "none") +
    labs(x = "Absolute delta PSI")

fig_c_title <- 
    ggdraw() + 
    draw_label(
	       "B cell stimulation leads to pervasive splicing effects (vs. Unstimulated cells)",
	       x = 0,
	       size = 9,
	       hjust = 0
	       ) +
    theme(text = element_text(size = 9),
	  plot.margin = margin(l = 1.25, unit = "lines"))

fig_c_grid <- 
    plot_grid(fig_c_title, volcano_plot, 
	      ncol = 1, rel_heights = c(.1, 1))


# Fig D
pathwaysgo <- gmtPathways("/lab-share/IM-Gutierrez-e2/Public/References/msigdb/c5.all.v2022.1.Hs.symbols.gmt.txt")
gobp <- keep(pathwaysgo, grepl("^GOBP", names(pathwaysgo)))

gsea_list <- 
    volcano_data |>
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

gsea_df <- 
    bind_rows(gsea_res, .id = "stim") |>
    as_tibble() |>
    group_by(stim) |>
    slice_max(n = 10, -log10(pval)) |>
    ungroup() |>
    mutate(stim = factor(stim, levels = c("TLR7c 24h", "BCRc 24h", "DN2c 24h")),
	   pathway = str_remove(pathway, "^GOBP_"),
	   pathway = str_replace(pathway, "POSITIVE_REGULATION", "POS_REGULATION"),
	   pathway = str_replace(pathway, "NEGATIVE_REGULATION", "NEG_REGULATION"),
	   pathway = str_trunc(pathway, 30))

gsea_plot <- 
    ggplot(gsea_df, 
	   aes(x = NES, 
	       y = reorder_within(pathway, by = NES, within = stim))) +
    geom_point(aes(fill = -log10(pval)), 
	       shape = 21, size = 2.5, stroke = .1) +
    scale_x_continuous(limits = c(1, 2),
		       breaks = c(1, 2)) +
    scale_y_reordered() +
    scale_fill_viridis_c(option = "magma", 
			 limits = c(1, 9),
			 breaks = c(1, 5, 9)) +
    facet_wrap(~stim, nrow = 1, scales = "free_y") +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 8),
	  axis.text.y = element_text(size = 6),
	  axis.title = element_text(size = 9),
	  strip.text = element_text(size = 9),
	  strip.clip = "off",
	  panel.grid.minor.x = element_blank(),
	  legend.position.inside = c(-.2, -.2),
	  legend.title = element_text(size = 9),
	  legend.text = element_text(size = 8)
	  ) +
    guides(fill = guide_colorbar(position = "inside", 
				 direction = "horizontal",
				 title.position = "top",
				 barheight = .25, barwidth = 6)) +
    labs(y = NULL)

fig_d_title <- 
    ggdraw() + 
    draw_label(
	       "Gene ontology biological processes enriched in differentially spliced genes",
	       x = 0,
	       size = 9,
	       hjust = 0
	       ) +
    theme(text = element_text(size = 9),
	  plot.margin = margin(l = 1.25, unit = "lines"))

fig_d_grid <- 
    plot_grid(fig_d_title, gsea_plot, 
	      ncol = 1, rel_heights = c(.1, 1))



# Fig E
load("../bcell_rnaseq/4-splicing/results/TLR7.Rdata")

tlr7_batf <- 
    counts[grepl("clu_17669_+", rownames(counts)), ] |>
    apply(2, function(x) x/sum(x)) |>
    as_tibble(rownames = "junction") |>
    pivot_longer(-junction) |>
    separate(name, c("stim", "donor", "replic"), sep = "_", remove = FALSE) |>
    group_by(stim, junction) |>
    summarise(value = mean(value)) |>
    ungroup() |>
    separate(junction, c("chr", "start", "end", "cluster"), sep = ":", convert = TRUE) |>
    select(stim, start, end, value)

rm(counts)

load("../bcell_rnaseq/4-splicing/results/BCR.Rdata")

bcr_batf <- 
    counts[grepl("clu_17669_+", rownames(counts)), ] |>
    apply(2, function(x) x/sum(x)) |>
    as_tibble(rownames = "junction") |>
    pivot_longer(-junction) |>
    separate(name, c("stim", "donor", "replic"), sep = "_", remove = FALSE) |>
    group_by(stim, junction) |>
    summarise(value = mean(value)) |>
    ungroup() |>
    separate(junction, c("chr", "start", "end", "cluster"), sep = ":", convert = TRUE) |>
    select(stim, start, end, value)

rm(counts)

load("../bcell_rnaseq/4-splicing/results/DN2.Rdata")

dn2_batf <- 
    counts[grepl("clu_17669_+", rownames(counts)), ] |>
    apply(2, function(x) x/sum(x)) |>
    as_tibble(rownames = "junction") |>
    pivot_longer(-junction) |>
    separate(name, c("stim", "donor", "replic"), sep = "_", remove = FALSE) |>
    group_by(stim, junction) |>
    summarise(value = mean(value)) |>
    ungroup() |>
    separate(junction, c("chr", "start", "end", "cluster"), sep = ":", convert = TRUE) |>
    select(stim, start, end, value)

batf_df <- bind_rows("TLR7" = tlr7_batf, "BCR" = bcr_batf, "DN2" = dn2_batf)















\rq
# Final figure
top_grid <-
    plot_grid(fig_a_grid, NULL, fig_b_grid, 
	      nrow = 1,
	      rel_widths = c(1, .1, .7),
	      labels = c("a", "", "b"), label_size = 12)

ggsave("./fig6.png", 
       plot_grid(top_grid, NULL, fig_c_grid, NULL, fig_d_grid, ncol = 1,
		 rel_heights = c(1, .05, 1, .05, 1),
		 labels = c("", "", "c", "", "d"), label_size = 12) +
       theme(plot.background = element_rect(color = "white", fill = "white")),
       width = 6.5, height = 6)








