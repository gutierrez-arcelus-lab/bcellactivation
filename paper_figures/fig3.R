# ==============================================================================
# Description: Generates Figure 3 (CITE-seq Single Cell Landscape).
#              Panel A: UMAP colored by stimulation condition.
#              Panel B: UMAP colored by Louvain cluster.
#              Panel C: UMAPs of selected ADT protein markers.
#              Panel D: Dot plot of RNA cluster markers.
#              Panel E: Bar plots of gated cell populations.
#              Panel F: UMAPs of disease-relevant RNA expression.
# ==============================================================================

library(tidyverse)
library(Seurat)
library(cowplot)
library(ggrepel)
library(scico)
library(ggridges)
library(ggrastr)

# functions
clip_expression <- function(values) { 
    
    clip_q <- quantile(values, 0.995)

    case_when(values >= clip_q ~ clip_q,
	      TRUE ~ values)
}

# Universal helper function to guarantee perfect title/label alignment
# Anchors text to the absolute top (y=1) and pushes it 15pt right to clear the letter label
create_title <- function(text) {
    ggdraw() + 
    draw_label(text, x = 0, y = 1, vjust = 1, hjust = 0, size = 7) + 
    theme(plot.margin = margin(t = 5, b = 5, l = 15, unit = "pt"))
}

# -----------------------------------------------------------------------------
# GLOBAL AESTHETICS & SETTINGS
# -----------------------------------------------------------------------------

# Set global ggplot theme to enforce max 7pt font size
theme_set(theme_minimal(base_size = 7))

stim_colors <- 
    read_tsv("./figure_colors_final.txt") |>
    unite("stim", c(Condition, Time), sep = " ") |>
    mutate(stim = paste0(stim, "h")) |>
    select(stim, Hex) |>
    deframe()

# meta data
features_df <- 
    "../04_citeseq/data/cellranger/1984/filtered_feature_bc_matrix/features.tsv.gz" |>
    read_tsv(col_names = c("gene_id", "gene_name", "phenotype"))

# data 
bcells <- 
    read_rds("../04_citeseq/2-processing/data/v4_seurat_qced.rds")

umap_df <- 
    read_tsv("../04_citeseq/2-processing/data/v4_umap_df.tsv") |>
    mutate(hto = factor(hto, levels = names(stim_colors)))

bcells@meta.data <- 
    bcells@meta.data |>
    mutate(dmm_hto_call = str_replace(dmm_hto_call, "IL4", "IL-4c"),
	   dmm_hto_call = str_replace(dmm_hto_call, "TLR7", "TLR7c"),
	   dmm_hto_call = str_replace(dmm_hto_call, "BCR", "BCRc"),
	   dmm_hto_call = str_replace(dmm_hto_call, "DN2", "DN2c"))

# -----------------------------------------------------------------------------
# Fig A: UMAP colored by condition
# -----------------------------------------------------------------------------
umap_stims <-
    ggplot(umap_df, aes(umap_1, umap_2)) +
    geom_point_rast(aes(fill = hto), 
		    size = .6, shape = 21, stroke = .05, color = "black",
		    raster.dpi = 600) +
    scale_fill_manual(values = stim_colors) +
    theme(axis.title = element_text(size = 7),
	  axis.text = element_text(size = 7),
	  axis.ticks = element_blank(),
	  panel.grid = element_blank(),
	  legend.text = element_text(size = 7, margin = margin(l = 0)),
	  legend.title = element_text(size = 7, margin = margin(b = 0)),
	  legend.margin = margin(0, 0, 0, 0),
	  legend.key.spacing.y = unit(-.2, "lines"),
	  plot.margin = margin(0, 0, 0, .5, unit = "lines")
	  ) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    guides(fill = guide_legend(title = "Stim:",
			       title.position = "top",
			       override.aes = list(size = 3)))

# -----------------------------------------------------------------------------
# Fig B: UMAP colored by cluster
# -----------------------------------------------------------------------------
cluster_labels <-
    umap_df |>
    mutate(cluster = factor(cluster)) |>
    group_by(cluster) |>
    summarise_at(vars(umap_1, umap_2), mean) |>
    ungroup()

cluster_colors <- c("#A8CDE2", "#3B83B9", "#E3362C", "#F9B56F", "#FC9230", "#DDA086",
		    "#9F7BB8", "#987898", "#F1E78D", "#B05D2F", "#83BF98", "#6ABD5D",
		    "#7E9D59", "#F4817F")

umap_clust <- 
    ggplot(umap_df, aes(umap_1, umap_2)) +
    geom_point_rast(aes(fill = factor(cluster)), 
		    size = .6, shape = 21, stroke = .05, color = "grey30",
		    raster.dpi = 600) +
    geom_label(data = cluster_labels |> filter(! cluster %in% c(12, 13)), 
	       aes(x = umap_1, y = umap_2, label = cluster),
	       label.padding = unit(0.1, "lines"),
	       size = 7, size.unit = "pt", alpha = .5, fontface = "bold") +
    geom_label_repel(data = cluster_labels |> filter(cluster %in% c(12, 13)), 
	       aes(x = umap_1, y = umap_2, label = cluster),
	       label.padding = unit(0.1, "lines"),
	       min.segment.length = 0, force_pull = 0, nudge_y = -1,
	       size = 7/.pt, alpha = 1, fontface = "bold") +
    scale_fill_manual(values = cluster_colors) +
    theme(
          axis.text = element_text(size = 7),
          axis.title = element_text(size = 7),
	  panel.grid = element_blank(),
	  plot.margin = margin(0, .5, 0, 0, unit = "lines")
	  ) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    guides(fill = "none") 

# -----------------------------------------------------------------------------
# Fig C: UMAPs of protein markers
# -----------------------------------------------------------------------------
marker_prots <- 
    c('IgM', 'IgD', 'CD1c', 
      'CD11c', 'CD21', 'CD23',
      'CD27', 'CD38', 'CD62L',
      'CD69', 'CD83', 'CD86')

marker_prots_df <-
    bcells@assays$ADT@data[marker_prots, ] |>
    as_tibble(rownames = "feature_id") |>
    pivot_longer(-feature_id, names_to = "barcode") |>
    left_join(select(umap_df, barcode, umap_1, umap_2), join_by(barcode)) |>
    mutate(feature_id = factor(feature_id, levels = marker_prots)) |>
    group_by(feature_id) |>
    mutate(value = clip_expression(value)) |>
    ungroup() 

umaps_prot_markers <-
    marker_prots_df |>
    group_split(feature_id) |>
    map(~ggplot(data = . |> arrange(value), aes(x = umap_1, y = umap_2)) +
	geom_point_rast(aes(color = value), 
			size = .1, stroke = 0,
			raster.dpi = 600) +
	scale_color_scico(palette = "lajolla",
			  label = function(x) sprintf("%.1f", x)) +
	facet_wrap(~feature_id, ncol = 1) +
	theme(axis.title = element_blank(),
	      axis.text = element_blank(),
	      legend.text = element_text(size = 7, margin = margin(l = 0)),
	      legend.margin = margin(0, 0, 0, -0.5, "lines"),
	      legend.box.spacing = unit(8, "pt"),
	      panel.grid = element_blank(),
	      strip.text = element_text(size = 7, margin = margin(t = 0, b = 0)),
	      strip.clip = "off",
	      panel.spacing = unit(0, "lines"),
	      plot.margin = margin(0, 0, 0, 0)
	      ) +
	labs(color = NULL) +
	guides(color = guide_colorbar(barwidth = .2, barheight = 2))
	) |>
    {function(x) plot_grid(plotlist = x, ncol = 3)}() +
    theme(plot.margin = margin(0, 0, 0, 0.5, unit = "lines"))

# -----------------------------------------------------------------------------
# Fig D: Dot plot of cluster markers
# -----------------------------------------------------------------------------
cluster_markers <- 
    read_tsv("../04_citeseq/2-processing/data/v4_cluster_markers.tsv")

top_markers <- 
    cluster_markers |>
    filter(p_val_adj <= 0.05, avg_log2FC > 0) |>
    group_by(cluster) |>
    top_n(10, abs(avg_log2FC)) |>
    ungroup()

top_markers_uniq <- 
    top_markers |>
    group_by(gene) |>
    slice_min(cluster) |>
    ungroup() |>
    arrange(cluster, p_val_adj, desc(abs(avg_log2FC))) |>
    select(gene, gene_name, gene_cluster = cluster)

curated_markers <-
    tribble(~cluster, ~gene_name,
	    0, "XBP1", 
	    0, "NR4A3", 
	    0, "IGHE", 
	    1, "IRF1",
	    1, "TBX21",
	    1, "GBP4",
	    2, "CCL3", 
	    2, "CCL4", 
	    3, "HILPDA",
	    4, "CD27",
	    4, "ITGAX",
	    5, "LMO2", 
	    5, "CD24", 
	    5, "FCER2", 
	    6, "SELL",
	    7, "TNFRSF13B",
	    8, "CCL17", 
	    9, "STMN1", 
	    10, "HMGB2",
	    10, "MKI67",
	    10, "CDC20",
	    11, "STAT1", 
	    11, "LGALS3", 
	    11, "IRF5",
	    12, "FOS", 
	    12, "DUSP1", 
	    12, "TCL1A", 
	    12, "TSC22D3",
	    13, "PRDM1", 
	    13, "IGHA1",
	    13, "JCHAIN",
	    13, "MZB1")

marker_genes <- 
    inner_join(curated_markers, cluster_markers) |>
    arrange(cluster, desc(avg_log2FC))

markers_plot_data <- 
    DotPlot(object = bcells, features = marker_genes$gene) |>
    {function(x) as_tibble(x$data)}() |>
    left_join(select(marker_genes, gene, gene_name), join_by(features.plot == gene)) |>
    select(cluster = id, gene = features.plot, gene_name, avg.exp, avg.exp.scaled, pct.exp) |>
    mutate(gene_name = factor(gene_name, levels = marker_genes$gene_name))

dotplot_markers <- 
    ggplot(markers_plot_data, aes(x = cluster, y = gene_name)) +
    geom_point(aes(size = pct.exp, fill = avg.exp.scaled), stroke = 0.2, shape = 21) +
    scale_fill_gradient2(low = "Light Sky Blue", mid = "lightyellow", high = "Dark Red",
			 midpoint = 0) +
    scale_size(range = c(0.1, 3), 
	       limits = c(0, 100),
	       breaks = c(0, 33, 66, 100)) +
    theme(
	  axis.text.x = element_text(size = 7),
	  axis.text.y = element_text(size = 7, face = 'italic'),
	  panel.grid = element_line(linewidth = .25, color = "grey90"),
	  legend.margin = margin(t = 0, r = 0, b = 0, l = -.25, unit = "lines"),
	  legend.title = element_text(size = 7),
	  legend.key.spacing.y = unit(-.5, "lines"),
	  plot.margin = margin(0, 0.5, 0, 0, unit = "lines")
	  ) +
    guides(fill = guide_colorbar(order = 1, position = "right",
				 barwidth = .25, barheight = 4)) +
    labs(x = NULL, y = NULL, fill = "Scaled\nExpression", size = "%\nExpressed")

# -----------------------------------------------------------------------------
# Fig E: Bar plots of cell populations
# -----------------------------------------------------------------------------
adt_df <- 
    bcells@assays$ADT@data[c('IgD', 'CD27', 'CD11c'), ] |>
    t() |>
    as_tibble(rownames = 'barcode') |>
    left_join(bcells@meta.data |> as_tibble(rownames = 'barcode') |> select(barcode, dmm_hto_call, seurat_clusters)) |>
    select(barcode, hto = dmm_hto_call, cluster = seurat_clusters, everything())

gates_genes <- 
    marker_genes |> 
    filter(gene_name %in% c('MKI67', 'XBP1')) |>
    select(gene, gene_name)

rna_df <-
    bcells@assays$RNA@data[gates_genes$gene, ] |>
    as.matrix() |>
    t() |>
    as_tibble(rownames = 'barcode') |>
    pivot_longer(-barcode, names_to = 'gene') |>
    left_join(gates_genes) |>
    select(-gene) |>
    pivot_wider(names_from = gene_name, values_from = value) |>
    left_join(bcells@meta.data |> as_tibble(rownames = 'barcode') |> select(barcode, dmm_hto_call, seurat_clusters)) |>
    select(barcode, hto = dmm_hto_call, cluster = seurat_clusters, everything())

subsets_tmp <- 
    left_join(adt_df, rna_df) |>
    mutate(
	   "IgD<sup>+</sup> CD27<sup>+</sup>"  = ifelse(IgD >= .75 & CD27 >= .25, 1, 0),
	   "IgD<sup>+</sup> CD27<sup>–</sup>"  = ifelse(IgD >= .75 & CD27 < .25, 1, 0),
	   "IgD<sup>–</sup> CD27<sup>+</sup>"  = ifelse(IgD < .75 & CD27 >= .25, 1, 0),
	   "IgD<sup>–</sup> CD27<sup>–</sup>"  = ifelse(IgD < .75 & CD27 < .25, 1, 0),
	   "CD11c<sup>+</sup>"                 = ifelse(CD11c >= .75, 1, 0),
	   "*XBP1* <sup>high</sup>"            = ifelse(XBP1 >= 2.5, 1, 0),
	   "*MKI67* <sup>high</sup>"           = ifelse(MKI67 >= 1, 1, 0))

subsets_df <- 
    subsets_tmp |>
    select(-(cluster:MKI67)) |>
    pivot_longer(-(barcode:hto), names_to = 'subtype') |>
    mutate(hto = factor(hto, names(stim_colors)),
           subtype = fct_inorder(subtype)) |>
    group_by(hto, subtype) |>
    summarise(pct = mean(value)) |>
    ungroup()

subsets_stim_donor <- 
    subsets_tmp |>
    select(-cluster, -IgD, -CD27, -CD11c, -XBP1, -MKI67) |>
    left_join(select(umap_df, barcode, donor_id), join_by(barcode)) |>
    select(barcode, donor_id, hto, everything()) |>
    pivot_longer(-(barcode:hto), names_to = "subtype") |>
    group_by(donor_id, hto, subtype) |>
    summarize(mean_prop = mean(value)) |>
    ungroup() |>
    mutate(hto = factor(hto, names(stim_colors)),
           subtype = factor(subtype, levels = levels(subsets_df$subtype)))

#subsets_stim_donor |>
#    group_by(hto, subtype) |>
#    summarise(min_prop = min(mean_prop),
#              max_prop = max(mean_prop)) |>
#    ungroup() 
#
#mean_min_max <- 
#    left_join(subsets_df, subsets_stim_donor) |>
#    mutate(subtype = str_remove_all(subtype, "\\*|<sup>|</sup>"),
#           hto = factor(hto, levels = levels(umap_df$hto))) |>
#    arrange(hto) |>
#    mutate_at(vars(pct, mean_prop), ~round(.x * 100, 2))
#
## In BCRc 72h, proportion of cells per cluster
#umap_df |>
#    filter(hto == "BCRc 72h") |>
#    count(cluster) |>
#    mutate(p = n/sum(n) * 100) |>
#    ungroup()
#
#umap_df |>
#    filter(hto == "BCRc 72h") |>
#    count(cluster, donor_id) |>
#    group_by(donor_id) |>
#    mutate(p = n/sum(n) * 100) |>
#    ungroup() |>
#    group_by(cluster) |>
#    summarise(min_p = min(p),
#              max_p = max(p)) |>
#    ungroup()
#
#umap_df |>
#    filter(hto == "BCRc 72h") |>
#    count(cluster, donor_id) |>
#    group_by(donor_id) |>
#    mutate(p = n/sum(n) * 100) |>
#    filter(cluster %in% 9:10) |>
#    summarise(p9_10 = sum(p)) |>
#    ungroup() |>
#    summarise(min_p = min(p9_10),
#              max_p = max(p9_10))

proportions_plot <- 
    ggplot(subsets_df, aes(x = hto, y = pct)) +
    geom_col(aes(color = hto, fill = hto), linewidth = .3, alpha = .25) +
    geom_jitter(data = subsets_stim_donor, 
                aes(x = hto, y = mean_prop, color = hto), 
                size = .6, width = .2) +
    scale_color_manual(values = stim_colors) +
    scale_fill_manual(values = stim_colors) +
    facet_wrap(~subtype, nrow = 1, scales = "free_y") +
    ggh4x::facetted_pos_scales(y = list(
					subtype == "IgD<sup>+</sup> CD27<sup>+</sup>" ~ scale_y_continuous(limits = c(0, 0.2), breaks = c(0, 0.2)),
					subtype == "IgD<sup>+</sup> CD27<sup>–</sup>" ~ scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)),
					subtype == "IgD<sup>–</sup> CD27<sup>+</sup>" ~ scale_y_continuous(limits = c(0, .6), breaks = c(0, 0.2, 0.4, 0.6)),
					subtype == "IgD<sup>–</sup> CD27<sup>–</sup>" ~ scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)),
					subtype == "CD11c<sup>+</sup>" ~ scale_y_continuous(limits = c(0, .3), breaks = c(0, 0.1, 0.2, 0.3)),
					subtype == "*XBP1* <sup>high</sup>" ~ scale_y_continuous(limits = c(0, .6), breaks = c(0, 0.2, 0.4, 0.6)),
					subtype == "*MKI67* <sup>high</sup>" ~ scale_y_continuous(limits = c(0, .2), breaks = c(0, 0.1, 0.2)))) +
    theme(
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 7),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.spacing.x = unit(.5, "lines"),
        axis.ticks.x = element_blank(),
        strip.text = ggtext::element_markdown(size = 7),
        strip.clip = "off",
        plot.margin = margin(t = 0.5, r = 0.5, b = 0, l = 0.5, unit = "lines")
	) +
    guides(color = "none", fill = "none") +
    labs(x = NULL, y = NULL)

# -----------------------------------------------------------------------------
# Fig F: UMAPs of disease-relevant genes
# -----------------------------------------------------------------------------
disease_genes <- 
    features_df |>
    filter(gene_name %in% c("IRF5", "IRF8", "MYC", "CXCR5", "CTSH", "RGS1", "ZBTB38"), 
	   phenotype == "Gene Expression") |>
    select(gene_name, gene_id)

dis_genes_expression <-
    bcells@assays$RNA@data[disease_genes$gene_id, ] |>
    as_tibble(rownames = "gene_id") |>
    left_join(disease_genes, join_by(gene_id)) |>
    pivot_longer(-c(gene_id, gene_name), names_to = "barcode") |>
    left_join(select(umap_df, barcode, umap_1, umap_2), join_by(barcode))

umaps_disease <-
    dis_genes_expression |>
    group_split(gene_name) |>
    map(function(x) 
	ggplot(data = x |> arrange(value) |> mutate(barcode = fct_inorder(barcode)), 
		aes(x = umap_1, y = umap_2)) +
	geom_point_rast(aes(color = value), 
			size = .3, stroke = 0,
			raster.dpi = 600) +
	scale_color_gradientn(colors = c("grey85","#FFF7EC","#FEE8C8","#FDD49E","#FDBB84",
					 "#FC8D59","#EF6548","#D7301F","#B30000","#7F0000"),
			      label = function(x) sprintf("%.1f", x)) +
	facet_wrap(~gene_name) +
	theme(axis.title = element_blank(),
	      axis.text = element_blank(),
	      legend.text = element_text(size = 7, margin = margin(l = 0)),
	      legend.margin = margin(0, 0, 0, -0.5, "lines"),
	      legend.box.spacing = unit(8, "pt"),
	      panel.grid = element_blank(),
	      strip.text = element_text(size = 7, face = 'italic', margin = margin(0, 0, 0, 0)),
	      strip.clip = "off",
	      panel.spacing = unit(0, "lines"),
	      plot.margin = margin(t = 0, r = 0, b = 0, l = 0)
	      ) +
	labs(color = NULL) +
	guides(color = guide_colorbar(barwidth = .2, barheight = 2))
	) |>
    {function(x) plot_grid(plotlist = x, nrow = 1)}() +
    theme(plot.margin = margin(0.5, .5, 0, 0.5, "lines"))

# -----------------------------------------------------------------------------
# Build final figure and export
# -----------------------------------------------------------------------------
fig_a_title <- create_title("UMAP shows separation of stimuli and time points")
fig_b_title <- create_title("Louvain clustering reveals B cell states")
fig_c_title <- create_title("Selected proteins markers")
fig_d_title <- create_title("RNA expression of marker genes for clusters in (b)")
fig_e_title <- create_title("Proportion of cell subsets in each condition in (a)")
fig_f_title <- create_title("RNA expression of selected disease risk genes in different B cell states")

# Row 1 (A & B)
guide_a <- get_plot_component(umap_stims, 'guide-box-right', return_all = TRUE)
plot_a_with_guide <- plot_grid(umap_stims + guides(fill = "none"), NULL, guide_a, nrow = 1, rel_widths = c(1, 0.05, 0.2))

titles_ab <- plot_grid(fig_a_title, NULL, fig_b_title, nrow = 1, rel_widths = c(1, 0.1, 0.8), 
                       labels = c('a', ' ', 'b'), label_size = 10, label_y = .9, vjust = 1)

plots_ab  <- plot_grid(plot_a_with_guide, NULL, umap_clust, nrow = 1, rel_widths = c(1, 0.1, 0.8))

grid_ab   <- plot_grid(titles_ab, plots_ab, ncol = 1, rel_heights = c(0.1, 1))

# Row 2 (C & D)
titles_cd <- plot_grid(fig_c_title, NULL, fig_d_title, nrow = 1, rel_widths = c(0.75, 0.05, 1), 
                       labels = c('c', ' ', 'd'), label_size = 10, label_y = .9, vjust = 1)

plots_cd  <- plot_grid(umaps_prot_markers, NULL, dotplot_markers, nrow = 1, rel_widths = c(0.75, 0.05, 1))

grid_cd   <- plot_grid(titles_cd, plots_cd, ncol = 1, rel_heights = c(0.1, 1))

# Row 3 (E)
title_e <- plot_grid(fig_e_title, labels = 'e', label_size = 10, label_y = .9, vjust = 1)
grid_e  <- plot_grid(title_e, proportions_plot, ncol = 1, rel_heights = c(0.1, 1))

# Row 4 (F)
title_f <- plot_grid(fig_f_title, labels = 'f', label_size = 10, label_y = .9, vjust = 1)
grid_f  <- plot_grid(title_f, umaps_disease, ncol = 1, rel_heights = c(0.1, 1))

# Final Combination
final_figure <- 
    plot_grid(grid_ab, NULL, grid_cd, NULL, grid_e, NULL, grid_f, 
          ncol = 1, rel_heights = c(0.75, 0.025, 1.05, 0.05, 0.3, 0.05, 0.33)) +
    theme(plot.background = element_rect(fill = "white", color = "white")) 

ggsave("./pdf/fig3.pdf", 
       final_figure,
       width = 179,
       height = 216,
       units = "mm",
       dpi = 600,
       device = cairo_pdf
       )
