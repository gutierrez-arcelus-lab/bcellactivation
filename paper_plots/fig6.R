library(tidyverse)
library(glue)
library(ggrepel)
library(tidytext)
library(tximport)
library(cowplot)
library(fgsea)
library(locuszoomr)

fig_colors <- 
    read_tsv("./figure_colors.txt", col_names = c("stim", "timep", "color")) |>
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
    read_tsv("../bcell_rnaseq/4-splicing/data/metadata_qced.tsv") |>
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
    geom_point(aes(fill = stim), size = 2.5, shape = 21, stroke = .2) +
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
read_leaf <- function(contrast) {
    
    sig <- 
        glue("../bcell_rnaseq/4-splicing/results/{contrast}_cluster_significance.txt") |>
        read_tsv() |>
	filter(status == "Success") |>
        separate(cluster, c("chr", "cluster"), sep = ":") |>
        select(cluster, p, padj = p.adjust, genes)

    eff <- 
        glue("../bcell_rnaseq/4-splicing/results/{contrast}_effect_sizes.txt") |>
        read_tsv() |>
        mutate(cluster = sub("^.*(clu.*)$", "\\1", intron)) |>
        select(cluster, logef, deltapsi)

    inner_join(sig, eff, join_by(cluster))
}

leaf_df <- 
    read_lines("../bcell_rnaseq/4-splicing/data/contrasts.txt") |>
    {function(x) setNames(x, x)}() |>
    map_dfr(read_leaf, .id = "contrast") |>
    mutate(contrast = str_replace(contrast, "unstday0", "Unstim 0h"),
	   contrast = str_replace(contrast, "TLR7", "TLR7c 24h"),
	   contrast = str_replace(contrast, "BCR", "BCRc 24h"),
	   contrast = str_replace(contrast, "DN2", "DN2c 24h"),
	   contrast = str_replace(contrast, "vs.", " vs. "))

leaf_summary <- 
    leaf_df |>
    filter(padj <= 0.05, !is.na(genes)) |>
    distinct(contrast, genes) |>
    count(contrast) |>
    separate(contrast, c("stim_a", "stim_b"), sep = " vs. ") |>
    mutate_at(vars(stim_a:stim_b), ~factor(., levels = names(rna_colors))) |>
    arrange(stim_a, stim_b)
    
leaf_summ_plot <- 
    ggplot(leaf_summary, aes(x = stim_a, y = stim_b)) +
    geom_tile(aes(fill = n)) +
    geom_text(aes(label = scales::comma(n)),
	      color = "black", fontface = "bold", size = 8, size.unit = "pt") +
    scale_x_discrete(position = "top") +
    scale_fill_gradient(low = "beige", high = "Dark Goldenrod") +
    theme_minimal() +
    theme(
	  panel.grid = element_blank(),
	  axis.text = element_text(size = 8),
	  axis.title = element_blank(),
	  axis.ticks.y = element_blank(),
	  legend.position = "none",
    )

fig_b_title <- 
    ggdraw() + 
    draw_label(
	       "Number of differentially spliced genes\n ",
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
    filter(grepl("^Unstim", contrast)) |>
    mutate(absd = abs(deltapsi)) |>
    group_by(contrast, genes) |>
    slice_min(p) |>
    slice_max(absd) |>
    ungroup() |>
    distinct(contrast, genes, p, absd) |>
    mutate(stim = sub("^Unstim 0h vs. (.+)$", "\\1", contrast)) |>
    mutate(stim = factor(stim, levels = names(rna_colors)[-1])) |>
    arrange(stim) |>
    mutate(contrast = sub("^(Unstim 0h) vs. (.+)$", "\\2 vs. \\1", contrast),
	   contrast = fct_inorder(contrast))

volcano_labels <- 
    volcano_data |>
    group_by(contrast) |>
    top_n(10, -log10(p)) |>
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
    facet_wrap(~contrast, nrow = 1) +
    theme_minimal() +
    theme(
	  axis.text = element_text(size = 8),
	  axis.title = element_text(size = 9),
	  strip.text = element_text(size = 9),
	  panel.grid.minor.y = element_blank(),
	  panel.grid.major.y = element_blank()) +
    guides(color = "none") +
    labs(x = "Absolute \u0394PSI")

fig_c_title <- 
    ggdraw() + 
    draw_label(
	       "B cell stimulation leads to pervasive splicing effects",
	       x = 0,
	       size = 9,
	       hjust = 0
	       ) +
    theme(text = element_text(size = 9),
	  plot.margin = margin(l = 1.25, unit = "lines"))

fig_c_grid <- 
    plot_grid(fig_c_title, volcano_plot, ncol = 1, rel_heights = c(.1, 1)) |>
    plot_grid(labels = "c", label_size = 12)


#### Test
#bcr_vs_dn2 <- leaf_df |>
#    filter(contrast == "BCRc 24h vs. DN2c 24h") |>
#    mutate(absd = abs(deltapsi)) |>
#    group_by(contrast, genes) |>
#    slice_min(p) |>
#    slice_max(absd) |>
#    ungroup() |>
#    distinct(contrast, genes, p, absd)
#
#bcr_vs_dn2_list <- 
#    bcr_vs_dn2 |>
#    mutate(p = -log10(p)) |>
#    separate_rows(genes, sep = ",") |>
#    filter(genes %in% unique(unlist(gobp))) |>
#    select(SYMBOL = genes, stat = p) |>
#    arrange(desc(stat)) |>
#    deframe()
#
#bcr_vs_dn2_res <- 
#    fgsea(pathways = gobp, 
#	  stats = bcr_vs_dn2_list, 
#	  scoreType = "pos", 
#	  nproc = future::availableCores())
#
#bcr_vs_dn2_res |>
#    as_tibble() |>
#    arrange(padj)
####

# Fig D
pathwaysgo <- gmtPathways("/lab-share/IM-Gutierrez-e2/Public/References/msigdb/c5.all.v2022.1.Hs.symbols.gmt.txt")
gobp <- keep(pathwaysgo, grepl("^GOBP", names(pathwaysgo)))

gsea_list <- 
    volcano_data |>
    mutate(p = -log10(p)) |>
    separate_rows(genes, sep = ",") |>
    filter(genes %in% unique(unlist(gobp))) |>
    select(contrast, SYMBOL = genes, stat = p) |>
    arrange(contrast, desc(stat)) |>
    {function(x) split(x, x$contrast)}() |>
    map(~select(., -contrast)) |>
    map(deframe)

gsea_res <- 
    map(gsea_list, 
	~fgsea(pathways = gobp, stats = ., scoreType = "pos", nproc = future::availableCores()))

gsea_df <- 
    bind_rows(gsea_res, .id = "contrast") |>
    as_tibble() |>
    mutate(contrast = factor(contrast, levels = names(gsea_list))) |>
    group_by(contrast) |>
    slice_max(n = 10, -log10(pval)) |>
    ungroup() |>
    mutate(
	   pathway = str_remove(pathway, "^GOBP_"),
	   pathway = str_replace(pathway, "POSITIVE_REGULATION", "POS_REGULATION"),
	   pathway = str_replace(pathway, "NEGATIVE_REGULATION", "NEG_REGULATION"),
	   pathway = str_trunc(pathway, 30))

gsea_plot <- 
    ggplot(gsea_df, 
	   aes(x = NES, 
	       y = reorder_within(pathway, by = NES, within = contrast))) +
    geom_point(aes(fill = padj), 
	       shape = 21, size = 2.5, stroke = .1) +
    scale_x_continuous(limits = c(1, 2),
		       breaks = c(1, 2)) +
    scale_y_reordered() +
    scale_fill_stepsn(name = "P_adj:",
		      breaks = c(.01, .05, .1, .2),
		      colors = c("#f03b20", "#feb24c", "#ffeda0", "grey")) +
    facet_wrap(~contrast, nrow = 1, scales = "free_y") +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 8),
	  axis.text.y = element_text(size = 6),
	  axis.title = element_text(size = 9),
	  strip.text = element_text(size = 9, margin = margin(l = -90, unit = "pt")),
	  strip.clip = "off",
	  panel.grid.minor.x = element_blank(),
	  legend.title = element_text(size = 9),
	  legend.text = element_text(size = 8),
	  legend.box.spacing = unit(6, "pt")
	  ) +
    guides(fill = guide_colorbar(barheight = 6, barwidth = .25)) +
    labs(x = "Normalized enrichment score", y = NULL)

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
    plot_grid(fig_d_title, gsea_plot, ncol = 1, rel_heights = c(.1, 1)) |>
    plot_grid(labels = "d", label_size = 12)


# Fig E
make_leaf_plot <- function(contrast, cluster_lab, title_plot) {

    source("./plot_cluster.R")
    
    glue("../bcell_rnaseq/4-splicing/results/{contrast}.Rdata") |>
        load()

    plot_cluster(cluster_to_plot = cluster_lab, 
		 main_title = title_plot,
		 exons_table = exons_table,
		 meta = meta,
		 cluster_ids = cluster_ids,
		 counts = counts,
		 introns = introns,
		 snp_pos = NA
		 ) +
    theme(legend.position = "none",
	  plot.margin = margin(0, 0, 0, 0),
	  plot.title = element_text(size = 9, hjust = .5)) +
    ylab(NULL)
}

cd86_grid <- make_leaf_plot("unstday0vs.DN2", "clu_9774_+", c("DN2c 24h", "Unstim 0h"))
zbtb38_grid <- make_leaf_plot("BCRvs.DN2", "clu_9831_+", c("DN2c 24h", "BCRc 24h"))


ah <- AnnotationHub::AnnotationHub()
ens_data <- ah[["AH98047"]]

filter <- dplyr::filter
select <- dplyr::select

track_cd86 <- 
    locus(gene = "CD86",
	  flank = c(1e3, 1e3),
	  ens_db = ens_data)

track_zbtb38 <- 
    locus(gene = "ZBTB38",
	  flank = c(1e3, 1e3),
	  ens_db = ens_data)


track_plot_cd86 <- 
    gg_genetracks(track_cd86, cex.text = .7) + 
    geom_rect(xmin = 122054000/1e6, 
	      xmax = 122093000/1e6, 
	      ymin = -.1, ymax = .5,
	      fill = NA, color = "red") +
    theme_minimal() +
    theme(
	  axis.text = element_blank(),
	  axis.title = element_blank(),
	  axis.line.x = element_blank(),
	  axis.ticks.x = element_blank(),
	  panel.grid = element_blank()) +
    coord_cartesian(clip = "off")

track_plot_zbtb38 <- 
    gg_genetracks(track_zbtb38, cex.text = .7, filter_gene_name = "ZBTB38") + 
    geom_rect(xmin = 141386000/1e6, 
	      xmax = 141404000/1e6, 
	      ymin = -.1, ymax = .5,
	      fill = NA, color = "red") +
    theme_minimal() +
    theme(
	  axis.text = element_blank(),
	  axis.title = element_blank(),
	  axis.line.x = element_blank(),
	  axis.ticks.x = element_blank(),
	  panel.grid = element_blank()) +
    coord_cartesian(clip = "off")

fig_e_title <- 
    ggdraw() + 
    draw_label(
	       "Differential splicing in the CD86 gene",
	       x = 0,
	       size = 9,
	       hjust = 0
	       ) +
    theme(text = element_text(size = 9),
	  plot.margin = margin(l = 1.25, unit = "lines"))

fig_e_grid <- 
    plot_grid(fig_e_title, track_plot_cd86, cd86_grid, 
	      ncol = 1, rel_heights = c(.1, .2, 1)) |>
    plot_grid(labels = "e", label_size = 12)

fig_f_title <- 
    ggdraw() + 
    draw_label(
	       "Differential splicing in the ZBTB38 gene",
	       x = 0,
	       size = 9,
	       hjust = 0
	       ) +
    theme(text = element_text(size = 9),
	  plot.margin = margin(l = 1.25, unit = "lines"))

fig_f_grid <- 
    plot_grid(fig_f_title, track_plot_zbtb38, zbtb38_grid, 
	      ncol = 1, rel_heights = c(.1, .2, 1)) |>
    plot_grid(labels = "f", label_size = 12)


bottom_grid <- plot_grid(fig_e_grid, fig_f_grid, nrow = 1)


# Final figure
top_grid <-
    plot_grid(fig_a_grid, NULL, fig_b_grid, 
	      nrow = 1,
	      rel_widths = c(1, .2, 1),
	      labels = c("a", "", "b"), label_size = 12)

ggsave("./fig6.png", 
       plot_grid(top_grid, NULL, fig_c_grid, NULL, fig_d_grid, NULL, bottom_grid,
		 ncol = 1,
		 rel_heights = c(1, .05, 1, .05, 1, .1, 1.5)) +
       theme(plot.background = element_rect(color = "white", fill = "white")),
       width = 6.5, height = 8.5, dpi = 300)

