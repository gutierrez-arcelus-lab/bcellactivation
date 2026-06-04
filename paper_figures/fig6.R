# ==============================================================================
# Description: Generates Figure 6 (Standard-input RNA-seq and Splicing).
#              Panel A: PCA of standard-input RNA-seq data.
#              Panel B: Number of differentially spliced genes.
#              Panel C: Volcano plots of splicing effects (deltaPSI).
#              Panel D: Gene ontology enrichment of differentially spliced genes.
#              Panel E: Differential splicing in the BCL2A1 gene.
#              Panel F: Differential splicing in the STAT6 gene.
#
# REQUIRED EXTERNAL DATASETS TO REPRODUCE THIS FIGURE:
# 1. MSigDB GO Biological Processes:
#    Download 'c5.all.v2022.1.Hs.symbols.gmt.txt' from the GSEA/MSigDB website.
# 2. Tian et al. Splicing Data (Nature Genetics, 2024):
#    Download Supplementary Table 11 (41588_2024_2019_MOESM4_ESM.xlsx).
# ==============================================================================

library(tidyverse)
library(glue)
library(ggrepel)
library(tidytext)
library(tximport)
library(cowplot)
library(fgsea)
library(locuszoomr)
library(GenomicRanges)

conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::select)
conflicted::conflicts_prefer(dplyr::slice)

# -----------------------------------------------------------------------------
# GLOBAL AESTHETICS & SETTINGS
# -----------------------------------------------------------------------------

# Set global ggplot theme to enforce max 7pt font size
theme_set(theme_minimal(base_size = 7))

# Helper for single-line titles
create_title <- function(text) {
    ggdraw() + 
    draw_label(text, x = 0, y = 1, vjust = 1, hjust = 0, size = 7) + 
    theme(plot.margin = margin(t = 5, b = 7, l = 15, unit = "pt"))
}

stim_colors <- 
    read_tsv("./figure_colors_final.txt") |>
    mutate(stim = glue("{Condition} {Time}h")) |>
    select(stim, Hex) |>
    deframe()

all_stims <- c("Unstim 0h", "TLR7c 24h", "BCRc 24h", "DN2c 24h")

# -----------------------------------------------------------------------------
# Fig A: PCA of standard-input RNA-seq
# -----------------------------------------------------------------------------
tx_to_gene <- 
    "../02_rnaseq_stdinput/2-quantification/data/gencode.v41.primary_assembly.annotation.gtf" |>
    rtracklayer::import(feature.type = "transcript") |>
    as_tibble() |>
    select(tx_id = transcript_id, gene_id, gene_name)

meta <- 
    read_tsv("../02_rnaseq_stdinput/1-mapping/data/metadata.tsv") |>
    separate(sample_id, c("stim", "donor", "replic"), sep = "_", remove = FALSE) |>
    select(sample_id, stim, donor) |>
    mutate(stim = recode(stim, 
             "unstday0" = "Unstim 0h",
             "TLR7" = "TLR7c 24h",
             "BCR" = "BCRc 24h",
             "DN2" = "DN2c 24h"),
       stim = factor(stim, levels = all_stims))

salmon_files <- 
    glue("../02_rnaseq_stdinput/2-quantification/results/{meta$sample_id}/quant.sf") |>
    set_names(meta$sample_id)

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
    geom_point(aes(fill = stim), size = 2.5, shape = 21, stroke = .5) +
    scale_fill_manual(values = stim_colors) +
    theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.title = element_text(size = 7),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 7),
      panel.grid = element_blank()) +
    guides(fill = guide_legend(override.aes = list(size = 4))) +
    labs(fill = "Stim:")


# -----------------------------------------------------------------------------
# Fig B: Number of differentially spliced genes
# -----------------------------------------------------------------------------
read_leaf <- function(contrast) {
    
    sig <- 
        glue("../02_rnaseq_stdinput/3-splicing/results/{contrast}_cluster_significance.txt") |>
        read_tsv() |>
        filter(status == "Success") |>
        separate(cluster, c("chr", "cluster"), sep = ":") |>
        select(cluster, p, padj = p.adjust, genes)

    eff <- 
        glue("../02_rnaseq_stdinput/3-splicing/results/{contrast}_effect_sizes.txt") |>
        read_tsv() |>
        mutate(cluster = sub("^.*(clu.*)$", "\\1", intron)) |>
        select(intron, cluster, logef, deltapsi)

    inner_join(sig, eff, join_by(cluster))
}

leaf_df <- 
    read_lines("../02_rnaseq_stdinput/3-splicing/data/contrasts.txt") |>
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
    mutate_at(vars(stim_a:stim_b), ~factor(., levels = all_stims)) |>
    arrange(stim_a, stim_b)
    
leaf_summ_plot <- 
    ggplot(leaf_summary, aes(x = stim_a, y = stim_b)) +
    geom_tile(aes(fill = n)) +
    geom_text(aes(label = scales::comma(n)),
          color = "black", fontface = "bold", size = 7, size.unit = "pt") +
    scale_x_discrete(position = "top") +
    scale_fill_gradient(low = "beige", high = "Dark Goldenrod") +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(size = 7),
      axis.title = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "none"
    )


# -----------------------------------------------------------------------------
# Fig C: Splicing effects (deltaPSI)
# -----------------------------------------------------------------------------
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
    mutate(stim = factor(stim, levels = all_stims[-1])) |>
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
		    size = 7 / .pt,
		    fontface = "italic",
		    # 1. Soften the lines so they don't overpower the plot
		    segment.size = 0.2,
		    segment.color = "grey60",
		    min.segment.length = 0,
		    # 3. Create an invisible bumper around each label and point
		    box.padding = 0.4,
		    point.padding = 0.25,
		    # 4. Aggressively push overlapping labels apart
		    force = 2,
		    max.overlaps = Inf) +
    scale_x_continuous(breaks = c(0, .5, 1),
		       labels = c("0", "0.5", "1")) +
    scale_color_manual(values = stim_colors) +
    facet_wrap(~contrast, nrow = 1) +
    theme_minimal() +
    theme(
	  axis.text = element_text(size = 7),
	  axis.title = element_text(size = 7),
	  strip.text = element_text(size = 7),
	  panel.grid.minor.y = element_blank(),
	  panel.grid.major.y = element_blank()) +
    guides(color = "none") +
    labs(x = "Absolute \u0394PSI")


# -----------------------------------------------------------------------------
# Fig D: Gene ontology enrichment
# -----------------------------------------------------------------------------
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
       pathway = str_replace(pathway, "POSITIVE_REGULATION", "Pos. Regulation"),
       pathway = str_replace(pathway, "NEGATIVE_REGULATION", "Neg. Regulation"),
       pathway = str_replace_all(pathway, "_", " "),
       pathway = str_to_title(pathway),
       pathway = str_replace_all(pathway, "\\bRna\\b", "RNA"),
       pathway = str_replace_all(pathway, "\\bDna\\b", "DNA"),
       pathway = str_replace_all(pathway, "\\bAtp\\b", "ATP"),
       pathway = str_replace_all(pathway, "\\bMrna\\b", "mRNA"),
       pathway = str_replace_all(pathway, "\\bOf\\b", "of"),
       pathway = str_trunc(pathway, 36))

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
    theme(axis.text.x = element_text(size = 7),
	  axis.text.y = element_text(size = 7),
	  axis.title = element_text(size = 7),
	  strip.text = element_text(size = 7, margin = margin(b = 7, l = -90, unit = "pt")),
	  strip.clip = "off",
	  panel.grid.minor.x = element_blank(),
	  legend.title = element_text(size = 7),
	  legend.text = element_text(size = 7),
	  legend.box.spacing = unit(6, "pt"),
	  legend.margin = margin(0, 0, 0, 0),
	  legend.box.margin = margin(0, 0, 0, 0)
    ) +
    guides(fill = guide_colorbar(barheight = 6, barwidth = .25)) +
    coord_cartesian(clip = "off") +
    labs(x = "Normalized enrichment score", y = NULL)


# -----------------------------------------------------------------------------
# Fig E & F: Differential Splicing Examples
# -----------------------------------------------------------------------------
make_leaf_plot <- function(contrast, cluster_lab, title_plot) {

    source("./plot_cluster.R")
    
    glue("../02_rnaseq_stdinput/3-splicing/results/{contrast}.Rdata") |>
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
      plot.title = element_text(size = 7, hjust = .5)) +
    ylab(NULL)
}

ah <- AnnotationHub::AnnotationHub()
ens_data <- ah[["AH98047"]]

conflicted::conflicts_prefer(dplyr::select)

tian <- 
    "../external/data/tian_splicing/41588_2024_2019_MOESM4_ESM.xlsx" |>
    readxl::read_excel(sheet = "Supplementary Table 11", skip = 2) |>
    rowid_to_column()

tian_tmp <- 
    tian |>
    filter(grepl("_RA_|_GD|SLE|Asthma|_AD|_LYM", trait)) |>
    select(rowid, Intron) |>
    separate(Intron, c("chrom", "start", "end", "clu"), sep = ":", convert = TRUE)

bcell_tmp <- 
    leaf_df |>
    rowid_to_column() |>
    group_by(cluster) |>
    filter(any(abs(deltapsi) >= 0.1)) |>
    ungroup() |>
    select(rowid, intron) |>
    separate(intron, c("chrom", "start", "end", "clu"), sep = ":", convert = TRUE)

tian_gr <- 
    GRanges(seqnames = tian_tmp$chrom,
            ranges = IRanges(start = tian_tmp$start, end = tian_tmp$end),
            rowid = tian_tmp$rowid)

bcell_gr <- 
    GRanges(seqnames = bcell_tmp$chrom,
            ranges = IRanges(start = bcell_tmp$start, end = bcell_tmp$end),
            rowid = bcell_tmp$rowid)
    
overlaps <- findOverlaps(bcell_gr, tian_gr, type = "equal")
   
matched_indices <- queryHits(overlaps)

my_matches <- 
    bcell_gr[matched_indices, ] |>
    as.data.frame() |>
    as_tibble()

bcl2a1_grid <- make_leaf_plot("unstday0vs.DN2", "clu_1476_-", c("DN2c 24h", "Unstim 0h"))
stat6_grid <- make_leaf_plot("unstday0vs.DN2", "clu_13231_-", c("DN2c 24h", "Unstim 0h"))

track_bcl2a1 <- 
    locus(gene = "BCL2A1",
          flank = c(1e3, 1e3),
          ens_db = ens_data)

track_stat6 <- 
    locus(gene = "STAT6",
          flank = c(1e3, 1e3),
          ens_db = ens_data)

track_plot_bcl2a1 <- 
    gg_genetracks(track_bcl2a1, cex.text = 7/12) + 
    theme_minimal() +
    theme(
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank()) +
    coord_cartesian(clip = "off")

track_plot_stat6 <- 
    gg_genetracks(track_stat6, cex.text = 7/12, filter_gene_name = "STAT6") + 
    geom_rect(xmin = 57109000/1e6,
              xmax = 57112200/1e6,
              ymin = -.002, ymax = .5,
              fill = NA, color = "red") +
    theme_minimal() +
    theme(
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank()) +
    coord_cartesian(clip = "off")


# -----------------------------------------------------------------------------
# STANDARDIZED ASSEMBLY & COMPILATION
# -----------------------------------------------------------------------------
fig_a_title <- create_title("Standard-input RNA-seq on 16 donors:\nPCA on gene expression data separates conditions")
fig_b_title <- create_title("Number of differentially spliced genes")
fig_c_title <- create_title("B cell stimulation leads to pervasive splicing effects")
fig_d_title <- create_title("Gene ontology biological processes enriched in differentially spliced genes")
fig_e_title <- create_title("Differential splicing in the BCL2A1 gene")
fig_f_title <- create_title("Differential splicing in the STAT6 gene")

# 1. Indent internal plots by 15pt to match titles and clear the labels
pca_indented       <- pca_plot + theme(plot.margin = margin(0.5, 0, 0, 0.5, "lines"))
leaf_summ_indented <- leaf_summ_plot + theme(plot.margin = margin(0.5, 0.5, 0, 0, "lines"))
volcano_indented   <- volcano_plot + theme(plot.margin = margin(0.5, 0.5, 0, 0.5, "lines"))
gsea_indented      <- gsea_plot + theme(plot.margin = margin(0.5, 0.5, 0, 0.5, "lines"))

bcl2a1_grid_indented <- 
    plot_grid(track_plot_bcl2a1, bcl2a1_grid, 
	      ncol = 1, rel_heights = c(.2, 1)) + 
    theme(plot.margin = margin(0, 0, 0, 0.5, "lines"))

stat6_grid_indented <- 
    plot_grid(track_plot_stat6, stat6_grid + theme(plot.margin = margin(t = 0, r = 0, b = .25, l = 0, "lines")),
	      ncol = 1, rel_heights = c(.2, 1)) + 
    theme(plot.margin = margin(0, 0.5, 0, 0, "lines"))

# 2. Rebuild the labeled grids over the indented plots (Labels sit at x=0)
fig_a_grid <- 
    plot_grid(fig_a_title, pca_indented, 
	      ncol = 1, rel_heights = c(.1, 1), 
	      labels = "a", label_size = 10, label_y = 1, vjust = 1)

fig_b_grid <- 
    plot_grid(fig_b_title, leaf_summ_indented, 
	      ncol = 1, rel_heights = c(.1, 1), 
	      labels = "b", label_size = 10, label_y = 1, vjust = 1)

fig_c_grid <- 
    plot_grid(fig_c_title, volcano_indented, 
	      ncol = 1, rel_heights = c(.1, 1), 
	      labels = "c", label_size = 10, label_y = 1, vjust = 1)

fig_d_grid <- 
    plot_grid(fig_d_title, gsea_indented, 
	      ncol = 1, rel_heights = c(.1, 1), 
	      labels = "d", label_size = 10, label_y = 1, vjust = 1)

fig_e_grid <- 
    plot_grid(fig_e_title, bcl2a1_grid_indented, 
	      ncol = 1, rel_heights = c(.1, 1), 
	      labels = "e", label_size = 10, label_y = 1, vjust = 1)

fig_f_grid <- 
    plot_grid(fig_f_title, stat6_grid_indented, 
	      ncol = 1, rel_heights = c(.1, 1), 
	      labels = "f", label_size = 10, label_y = 1, vjust = 1)

# 3. Assemble the main rows
top_grid    <- plot_grid(fig_a_grid, NULL, fig_b_grid, nrow = 1, rel_widths = c(1, .1, 1))
bottom_grid <- plot_grid(fig_e_grid, fig_f_grid, nrow = 1)

# 4. Final Combination
final_figure <- 
    plot_grid(top_grid, NULL, fig_c_grid, NULL, fig_d_grid, NULL, bottom_grid,
          ncol = 1,
          rel_heights = c(1, .05, 1, .05, 1, .1, 1.5)) +
    theme(plot.background = element_rect(color = "white", fill = "white"),
          plot.margin = margin(10, 0, -15, 0, "pt"))

# Final fig (Exported as 179mm PDF)
ggsave("./pdf/fig6.pdf", 
       final_figure,
       width = 179, 
       height = 216, 
       units = "mm",
       dpi = 600,
       device = cairo_pdf)
