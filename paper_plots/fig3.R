library(tidyverse)
library(Seurat)
library(cowplot)

# functions
clip_expression <- function(values) { 
    
    clip_q <- quantile(values, 0.99)

    case_when(values >= clip_q ~ clip_q,
	      TRUE ~ values)
}


# colors
stim_colors <- 
    read_tsv("../figure_colors.txt", col_names = c("stim", "time", "color")) |>
    mutate(stim = recode(stim, 
			 "IL4" = "IL-4c", 
			 "TLR7" = "TLR7c",
			 "BCR" = "BCRc", 
			 "DN2" = "DN2c")) |>
    unite("stim", c(stim, time), sep = " ") |>
    mutate(stim = paste0(stim, "h")) |>
    deframe()

# meta data
cellranger_dir_1984 <- 
    file.path("/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/B_cells_citeseq",
	      "SN0268787/broad/hptmp/curtism/bwh10x/KW10456_Maria",
	      "221101_10X_KW10456_bcl/cellranger-6.1.1/GRCh38/BRI-1984/outs",
	      "filtered_feature_bc_matrix")

features_df <- file.path(cellranger_dir_1984, "features.tsv.gz") |>
    read_tsv(col_names = c("gene_id", "gene_name", "phenotype"))

# data 
bcells <- read_rds("../citeseq/data/seuratv4_qced.rds")

umap_df <- 
    read_tsv("../citeseq/data/v4_umap_df.tsv") |>
    mutate(hto = str_replace(hto, "IL4", "IL-4c"),
	   hto = str_replace(hto, "TLR7", "TLR7c"),
	   hto = str_replace(hto, "BCR", "BCRc"),
	   hto = str_replace(hto, "DN2", "DN2c")) |>
    mutate(hto = factor(hto, levels = names(stim_colors)))

bcells@meta.data <- bcells@meta.data |>
    mutate(dmm_hto_call = str_replace(dmm_hto_call, "IL4", "IL-4c"),
	   dmm_hto_call = str_replace(dmm_hto_call, "TLR7", "TLR7c"),
	   dmm_hto_call = str_replace(dmm_hto_call, "BCR", "BCRc"),
	   dmm_hto_call = str_replace(dmm_hto_call, "DN2", "DN2c"))


# Figure A ####################################################################
umap_stims <-
    ggplot(umap_df, aes(umap_1, umap_2)) +
    geom_point(aes(fill = hto), size = .6, 
	       shape = 21, stroke = .05, color = "black") +
    scale_fill_manual(values = stim_colors) +
    theme_minimal() +
    theme(axis.title = element_text(size = 8),
	  axis.text = element_text(size = 8),
	  axis.ticks = element_blank(),
	  panel.grid = element_blank(),
	  legend.text = element_text(size = 8, margin = margin(l = 0)),
	  legend.title = element_text(size = 9, margin = margin(b = 0)),
	  legend.margin = margin(0, 0, 0, 0),
	  legend.key.spacing.y = unit(-.2, "lines"),
	  plot.margin = margin(0, 0, 0, 0)
	  ) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    guides(fill = guide_legend(title = "Stim:",
			       title.position = "top",
			       override.aes = list(size = 3)))


# Figure B ####################################################################
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
    geom_point(aes(fill = factor(cluster)), size = .6, 
	       shape = 21, stroke = .05, color = "grey30") +
    geom_label(data = cluster_labels |> filter(! cluster %in% c(12, 13)), 
	       aes(x = umap_1, y = umap_2, label = cluster),
	       label.padding = unit(0.1, "lines"),
	       size = 8, size.unit = "pt", alpha = .5, fontface = "bold") +
    ggrepel::geom_label_repel(data = cluster_labels |> filter(cluster %in% c(12, 13)), 
	       aes(x = umap_1, y = umap_2, label = cluster),
	       label.padding = unit(0.1, "lines"),
	       min.segment.length = 0, force_pull = 0, nudge_y = -1,
	       size = 3, alpha = 1, fontface = "bold") +
    scale_fill_manual(values = cluster_colors) +
    theme_minimal() +
    theme(
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 8),
	  panel.grid = element_blank(),
	  plot.margin = margin(0, 0, 0, 0)
	  ) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    guides(fill = "none") 

# Make top panel 
fig_a_title <- 
    ggdraw() + 
    draw_label(
	       "UMAP shows separation of stimuli and time points",
	       x = 0,
	       size = 9,
	       hjust = 0
	       ) +
    theme(text = element_text(size = 9),
	  plot.margin = margin(l = 1.25, unit = "lines"))

fig_b_title <- 
    ggdraw() + 
    draw_label(
	       "Louvain clustering reveals B cell states",
	       x = 0,
	       size = 9,
	       hjust = 0
	       ) +
    theme(text = element_text(size = 9),
	  plot.margin = margin(l = 1.25, unit = "lines"))

top_umaps <- 
    plot_grid(umap_stims + guides(fill = "none"), 
	      NULL, 
	      get_plot_component(umap_stims, 'guide-box-right', return_all = TRUE),
	      NULL,
	      umap_clust,
	      nrow = 1, rel_widths = c(1, .025, .2, .2, 1))

top_titles <- 
    plot_grid(fig_a_title, fig_b_title, nrow = 1, rel_widths = c(1, .8),
	      labels = c('a', 'b'), label_size = 12)


grid_ab <- plot_grid(top_titles, top_umaps, ncol = 1, rel_heights = c(.1, 1))


# Figure C ####################################################################
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
    map(~ggplot(data = ., aes(x = umap_1, y = umap_2)) +
	geom_point(aes(color = value), size = .2, stroke = 0) +
	scale_color_viridis_c(option = "magma") +
	facet_wrap(~feature_id, ncol = 1) +
	theme_minimal() +
	theme(axis.title = element_blank(),
	      axis.text = element_blank(),
	      legend.text = element_blank(),
	      legend.margin = margin(0, 0, 0, -0.5, "lines"),
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
    theme(plot.margin = margin(0, 0, 0, 0))

fig_c_title <- 
    ggdraw() + 
    draw_label(
	       "Selected proteins markers",
	       x = 0,
	       size = 9,
	       hjust = 0
	       ) +
    theme(text = element_text(size = 9),
	  plot.margin = margin(l = 1.25, unit = "lines"))

# Figure D ###################################################################
cluster_markers <- read_tsv("../citeseq/data/v4_cluster_markers.tsv")

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
    theme_minimal() +
    theme(
	  axis.text.x = element_text(size = 8),
	  axis.text.y = element_text(size = 8, face = 'italic'),
	  panel.grid = element_line(linewidth = .25, color = "grey90"),
	  legend.margin = margin(t = 0, r = 0.2, b = 0, l = -.5, unit = "lines"),
	  legend.title = element_text(size = 8),
	  legend.key.spacing.y = unit(-.5, "lines"),
	  plot.margin = margin(0, 0.1, 0, 0, unit = "lines")
	  ) +
    guides(fill = guide_colorbar(order = 1, position = "right",
				 barwidth = .25, barheight = 4)) +
    labs(x = NULL, y = NULL, fill = "Scaled\nExpression", size = "%\nExpressed")


fig_d_title <- 
    ggdraw() + 
    draw_label(
	       "RNA expression of marker genes for clusters in (b)",
	       x = 0,
	       size = 9,
	       hjust = 0
	       ) +
    theme(text = element_text(size = 9),
	  plot.margin = margin(l = 1.25, unit = "lines"))

# Middle panel

titles_cd <- 
    plot_grid(fig_c_title, NULL, fig_d_title, nrow = 1, rel_widths = c(.66, .05, 1),
	      labels = c('c', '', 'd'), label_size = 12)

fig_cd <- plot_grid(umaps_prot_markers, NULL, dotplot_markers, nrow = 1, rel_widths = c(.7, .1, 1))

grid_cd <- plot_grid(titles_cd, fig_cd, ncol = 1, rel_heights = c(.1, 1))


# Figure E ####################################################################
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

#testp <-
#    ggplot(rna_df, aes(x = cluster, y = MKI67)) +
#    geom_violin() +
#    geom_hline(yintercept = .0)
#
#
#testp <-
#    ggplot(adt_df, aes(x = cluster, y = IgD)) +
#    geom_violin() +
#    geom_hline(yintercept = .75)
#
#testp <-
#    ggplot(adt_df, aes(x = cluster, y = CD27)) +
#    geom_violin() +
#    geom_hline(yintercept = .25)
#
#testp <-
#    ggplot(adt_df, aes(x = cluster, y = CD11c)) +
#    geom_violin() +
#    geom_hline(yintercept = .25)
#
#ggsave("./testp.png", testp, height = 4)
#
subsets_df <- 
    left_join(adt_df, rna_df) |>
    mutate(
	   "IgD<sup>+</sup>  CD27<sup>+</sup>"  = ifelse(IgD >= .75 & CD27 >= .75, 1, 0),
	   "IgD<sup>–</sup>  CD27<sup>–</sup>" = ifelse(IgD < .75 & CD27 < .75, 1, 0),
	   "IgD<sup>–</sup>  CD27<sup>+</sup>" = ifelse(IgD < .75 & CD27 >= .75, 1, 0),
	   "IgD<sup>+</sup>  CD27<sup>–</sup>" = ifelse(IgD >= .75 & CD27 < .75, 1, 0),
	   "CD11c<sup>+</sup>" = ifelse(CD11c >= .75, 1, 0),
	   "*XBP1*<sup>high</sup>" = ifelse(XBP1 >= 2.5, 1, 0),
	   "*MKI67*<sup>+</sup>" = ifelse(MKI67 >= 1, 1, 0)) |>
    select(-(cluster:MKI67)) |>
    pivot_longer(-(barcode:hto), names_to = 'subtype') |>
    mutate(hto = factor(hto, names(stim_colors)),
	   subtype = fct_inorder(subtype)) |>
    group_by(hto, subtype) |>
    summarise(pct = mean(value)) |>
    ungroup()

#bcells_v5 <- read_rds("../citeseq/data/seurat_qced.rds")
#
#gating_subtypes <-
#    bcells_v5@meta.data |>
#    as_tibble(rownames = "barcode") |>
#    select(barcode, hto = dmm_hto_call, starts_with("is.pure_")) |>
#    pivot_longer(-(barcode:hto), names_to = "subtype") |>
#    mutate(subtype = str_remove(subtype, "is.pure_"),
#	   value = recode(value, "Impure" = "No", "Pure" = "Yes"))
#
#gating_scores <- 
#    bcells_v5@meta.data |>
#    as_tibble(rownames = "barcode") |>
#    select(barcode, hto = dmm_hto_call, ends_with("UCell")) |>
#    pivot_longer(-(barcode:hto), names_to = "subtype", values_to = "score") |>
#    mutate(subtype = str_remove(subtype, "_UCell"),
#	   subtype = recode(subtype, "memory" = "mem", "prolif." = "prolif"))
#
#gating_df <-
#    left_join(gating_subtypes, gating_scores, join_by(barcode, hto, subtype)) 
#
#gating_df <-
#    left_join(gating_subtypes, gating_scores, join_by(barcode, hto, subtype)) |>
#    group_by(subtype) |>
#    mutate(value = case_when(subtype %in% c("naive", "plasma") & value == "Yes" & score >= quantile(score, .8) ~ "Yes",
#			     subtype %in% c("naive", "plasma") & value == "Yes" & score < quantile(score, .8) ~ "No",
#			     .default = value)) |>
#    ungroup()
#
#subsets_df <-
#    gating_df |>
#    select(barcode, subtype, hto, value) |>
#    mutate(hto = factor(hto, levels = names(stim_colors))) |>
#    group_by(hto, subtype) |>
#    summarise(pct = mean(value == "Yes")) |>
#    ungroup() |>
#    mutate(subtype = recode(subtype, 
#			    "abc" = "CD11c<sup>+</sup>",
#			    "mem" = "IgD<sup>–</sup>  CD27<sup>+</sup>",
#			    "naive" = "IgD<sup>+</sup>  CD27<sup>–</sup>",
#			    "plasma" = "*XBP1*<sup>high</sup>",
#			    "prolif" = "*MKI67*<sup>+</sup>")) |>
#    arrange(hto, desc(pct)) |>
#    mutate(subtype = fct_inorder(subtype))
#
proportions_plot <- 
    ggplot(subsets_df, aes(x = hto, y = pct)) +
    geom_col(aes(fill = hto)) +
    scale_fill_manual(values = stim_colors) +
    facet_wrap(~subtype, nrow = 1, scales = "free_y") +
    ggh4x::facetted_pos_scales(y = list(
					subtype == "IgD<sup>+</sup>  CD27<sup>–</sup>" ~ scale_y_continuous(limits = c(0, 0.8), breaks = c(0, 0.4, 0.8)),
					subtype == "IgD<sup>–</sup>  CD27<sup>+</sup>" ~ scale_y_continuous(limits = c(0, .2), breaks = c(0, 0.1, 0.2)),
					subtype == "IgD<sup>+</sup>  CD27<sup>+</sup>" ~ scale_y_continuous(limits = c(0, 0.005), breaks = c(0, 0.005)),
					subtype == "IgD<sup>–</sup>  CD27<sup>–</sup>" ~ scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)),
					subtype == "CD11c<sup>+</sup>" ~ scale_y_continuous(limits = c(0, .2), breaks = c(0, 0.1, 0.2)),
					subtype == "*XBP1*<sup>high</sup>" ~ scale_y_continuous(limits = c(0, .5), breaks = c(0, 0.25, 0.5)),
					subtype == "*MKI67*<sup>+</sup>" ~ scale_y_continuous(limits = c(0, .2), breaks = c(0, 0.1, 0.2)))) +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
	  panel.grid.major.x = element_blank(),
	  panel.grid.minor.y = element_blank(),
	  panel.spacing.x = unit(1, "lines"),
	  axis.ticks.x = element_blank(),
	  strip.text = ggtext::element_markdown(size = 8),
	  strip.clip = "off",
	  plot.margin = margin(t = 5.5, r = 5.5/2, b = 5.5, l = 5.5/2, "pt")
	  ) +
    guides(fill = "none") +
    labs(x = NULL, y = NULL)


#proportions_plot <- 
#    ggplot(subsets_df, aes(x = hto, y = pct)) +
#    geom_col(aes(fill = hto)) +
#    scale_fill_manual(values = stim_colors) +
#    facet_wrap(~subtype, nrow = 1, scales = "free_y") +
#    ggh4x::facetted_pos_scales(y = list(subtype == "IgD<sup>+</sup>  CD27<sup>–</sup>" ~ scale_y_continuous(breaks = c(0, 0.4, 0.8)),
#					subtype == "IgD<sup>–</sup>  CD27<sup>+</sup>" ~ scale_y_continuous(breaks = c(0, 0.1, 0.2)),
#					subtype == "CD11c<sup>+</sup>" ~ scale_y_continuous(limits = c(0, .2), breaks = c(0, 0.1, 0.2)),
#					subtype == "*XBP1*<sup>high</sup>" ~ scale_y_continuous(breaks = c(0, 0.25, 0.5)),
#					subtype == "*MKI67*<sup>+</sup>" ~ scale_y_continuous(limits = c(0, .2), breaks = c(0, 0.1, 0.2)))) +
#    theme_minimal() +
#    theme(axis.text.x = element_blank(),
#	  panel.grid.major.x = element_blank(),
#	  panel.grid.minor.y = element_blank(),
#	  panel.spacing.x = unit(1, "lines"),
#	  axis.ticks.x = element_blank(),
#	  strip.text = ggtext::element_markdown(size = 8),
#	  strip.clip = "off",
#	  plot.margin = margin(t = 5.5, r = 5.5/2, b = 5.5, l = 5.5/2, "pt")
#	  ) +
#    guides(fill = "none") +
#    labs(x = NULL, y = NULL)


fig_e_title <- 
    ggdraw() + 
    draw_label(
	       "Proportion of cell subsets in each condition in (a)",
	       x = 0,
	       size = 9,
	       hjust = 0
	       ) +
    theme(text = element_text(size = 9),
	  plot.margin = margin(l = 1.25, unit = "lines"))

fig_e <- plot_grid(fig_e_title,
		   proportions_plot, 
		   ncol = 1, 
		   rel_heights = c(.1, 1),
		   labels = 'e', label_size = 12, vjust = .75)

# Figure F ####################################################################
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
    left_join(select(umap_df, barcode, umap_1, umap_2), join_by(barcode)) #|>
    #group_by(gene_id, gene_name) |>
    #mutate(value = clip_expression(value)) |>
    #ungroup()

umaps_disease <-
    dis_genes_expression |>
    group_split(gene_name) |>
    map(function(x) 
	ggplot(data = x |> arrange(value) |> mutate(barcode = fct_inorder(barcode)), 
		aes(x = umap_1, y = umap_2)) +
	geom_point(aes(color = value), size = .3, stroke = 0) +
	scale_color_gradientn(colors = c("grey85","#FFF7EC","#FEE8C8","#FDD49E","#FDBB84",
					 "#FC8D59","#EF6548","#D7301F","#B30000","#7F0000")) +
	facet_wrap(~gene_name) +
	theme_minimal() +
	theme(axis.title = element_blank(),
	      axis.text = element_blank(),
	      legend.text = element_blank(),
	      legend.margin = margin(0, 0, 0, -0.5, "lines"),
	      panel.grid = element_blank(),
	      strip.text = element_text(size = 8, face = 'italic', margin = margin(0, 0, 0, 0)),
	      strip.clip = "off",
	      panel.spacing = unit(0, "lines"),
	      plot.margin = margin(0, 0, 0, 0)
	      ) +
	labs(color = NULL) +
	guides(color = guide_colorbar(barwidth = .2, barheight = 2))
	) |>
    {function(x) plot_grid(plotlist = x, nrow = 1)}() +
    theme(plot.margin = margin(5.5, 0, 5.5, 0, "pt"))

fig_f_title <- 
    ggdraw() + 
    draw_label(
	       "RNA expression of selected disease risk genes in different B cell states",
	       x = 0,
	       size = 9,
	       hjust = 0
	       ) +
    theme(text = element_text(size = 9),
	  plot.margin = margin(l = 1.25, unit = "lines"))

fig_f <- plot_grid(fig_f_title, umaps_disease, ncol = 1, rel_heights = c(.1, 1),
		   labels = 'f', label_size = 12, vjust = 1)


ggsave("./fig3.png", 
       plot_grid(grid_ab, NULL, grid_cd, NULL, fig_e, NULL, fig_f, 
		 ncol = 1, rel_heights = c(.75, .025, 1.05, 0.05, .3, 0.05, .36)), 
       width = 6.5, height = 8.5, dpi = 600)

