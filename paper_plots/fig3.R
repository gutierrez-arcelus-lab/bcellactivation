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
    read_tsv("../figure_colors.txt", col_names = c("stim", "color")) |>
    mutate(stim = str_replace(stim, "_", " "),
	   stim = paste0(stim, "h")) |>
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
bcells <- read_rds("../citeseq/data/seurat_qced.rds")

umap_df <- 
    read_tsv("../citeseq/data/umap_df.tsv") |>
    mutate(hto = factor(hto, levels = names(stim_colors)))

# Figure A ####################################################################
umap_stims <-
    ggplot(umap_df, aes(umap_1, umap_2)) +
    geom_point(aes(fill = hto), size = .6, 
	       shape = 21, stroke = .05, color = "grey30") +
    scale_fill_manual(values = stim_colors) +
    theme_minimal() +
    theme(axis.title = element_blank(),
	  axis.text = element_blank(),
	  axis.ticks = element_blank(),
	  panel.grid = element_blank(),
	  legend.text = element_text(size = 8, margin = margin(l = 0)),
	  legend.title = element_text(size = 9, margin = margin(b = 0)),
	  legend.margin = margin(l = -1, unit = "lines"),
	  legend.key.spacing.y = unit(-.2, "lines"),
	  plot.margin = margin(r = 0)
	  ) +
    guides(fill = guide_legend(title = "Stim:",
			       title.position = "top",
			       override.aes = list(size = 3)))

fig_a_title <- 
    ggdraw() + 
    draw_label(
	       "UMAP shows separation of stimuli and\ntime points",
	       x = 0,
	       size = 9,
	       hjust = 0
	       ) +
    theme(text = element_text(size = 9),
	  plot.margin = margin(l = 1.25, unit = "lines"))

fig_a <- plot_grid(fig_a_title, umap_stims, ncol = 1, rel_heights = c(.1, 1))


# Figure B ####################################################################
marker_prots <- 
    c("IgM", "IgD", "CD62L", "CD38",
      "CD1c", "CD27", "CD11c", "CD86",
      "CD21", "CD23", "CD69", "CD83")

marker_prots_df <- 
    LayerData(bcells, assay = "ADT", layer = "data")[marker_prots, ] |>
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
	geom_point(aes(color = value), size = .1, stroke = 0) +
	scale_color_gradient(low = "grey90", high = "darkred") +
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
    {function(x) plot_grid(plotlist = x, ncol = 4)}() +
    theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))

fig_b_title <- 
    ggdraw() + 
    draw_label(
	       "Surface proteins markers",
	       x = 0,
	       size = 9,
	       hjust = 0
	       ) +
    theme(text = element_text(size = 9),
	  plot.margin = margin(l = 1.25, unit = "lines"))

fig_b <- plot_grid(fig_b_title, umaps_prot_markers, ncol = 1, rel_heights = c(.1, 1))


# Figure C ####################################################################
cluster_labels <-
    umap_df |>
    mutate(cluster = factor(cluster)) |>
    group_by(cluster) |>
    summarise_at(vars(umap_1, umap_2), mean) |>
    ungroup()

umap_clust <- 
    ggplot(umap_df, aes(umap_1, umap_2)) +
    geom_point(aes(fill = factor(cluster)), size = .6, 
	       shape = 21, stroke = .05, color = "grey30") +
    geom_label(data = cluster_labels, 
	       aes(x = umap_1, y = umap_2, label = cluster),
	       label.padding = unit(0.1, "lines"),
	       size = 8, size.unit = "pt", alpha = .5, fontface = "bold") +
    scale_fill_manual(values = pals::kelly(n = length(unique(umap_df$cluster)) + 1 )[-1] ) +
    theme_minimal() +
    theme(
          axis.text = element_blank(),
          axis.title = element_blank(),
	  panel.grid = element_blank(),
	  plot.margin = margin(t = 5.5, r = 0, b = 5.5, l = 0, unit = "pt")
	  ) +
    guides(fill = "none") 

fig_c_title <- 
    ggdraw() + 
    draw_label(
	       "Louvain clustering",
	       x = 0,
	       size = 9,
	       hjust = 0
	       ) +
    theme(text = element_text(size = 9),
	  plot.margin = margin(l = 1.25, unit = "lines"))

fig_c <- plot_grid(fig_c_title, umap_clust, ncol = 1, rel_heights = c(.1, 1))







# Figure D ###################################################################

cluster_markers <- read_tsv("../citeseq/data/cluster_markers.tsv")

top_markers <- 
    cluster_markers |>
    filter(p_val_adj <= 0.05) |>
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
    c("XBP1", "JCHAIN", "IGHE", "FCER2", "NR4A3",
      "IRF1", "TBX21", "GBP4", "STAT1", "SLAMF7",
      "CCL4",
      "IGHA1", "ITGAX", "CD27", "TCL1A", "CD24", "TNFRSF13B", "CD1C", 
      "MKI67", "CDC20", "UHRF1", 
      "FOS", "DUSP1")

marker_genes <- 
    cluster_markers |>
    filter(gene_name %in% curated_markers) |>
    group_by(gene_name) |>
    slice_min(cluster) |>
    ungroup() |>
    arrange(cluster, p_val_adj, desc(abs(avg_log2FC)))

markers_expression <- 
    LayerData(bcells, assay = "RNA", layer = "data")[marker_genes$gene, ] |>
    as_tibble(rownames = "gene") |>
    left_join(select(marker_genes, gene, gene_name, cluster), join_by(gene)) |>
    mutate(gene_name = factor(gene_name, levels = curated_markers)) |>
    pivot_longer(-c(gene, gene_name, cluster), names_to = "barcode") |>
    rename("gene_cluster" = cluster) |>
    left_join(as_tibble(bcells@meta.data, rownames = "barcode") |> select(barcode, seurat_clusters),
	      join_by(barcode)) |>
    select(gene_name, gene_cluster, barcode, seurat_clusters, value)

markers_summary <- 
    markers_expression |>
    group_by(gene_name, gene_cluster, seurat_clusters) |>
    summarise(pct = mean(value > 0),
	      mean_value = mean(value[value > 0])) |>
    ungroup() |>
    mutate(mean_value = replace_na(mean_value, 0)) |>
    group_by(gene_name) |>
    mutate(mean_value = mean_value/max(mean_value)) |>
    ungroup() |>
    mutate(seurat_clusters = paste0("C-", seurat_clusters),
	   seurat_clusters = factor(seurat_clusters, levels = paste0("C-", 0:13)))


dotplot_markers <- 
    ggplot(markers_summary, aes(x = gene_name, y = seurat_clusters)) +
    geom_point(aes(size = pct, fill = mean_value), stroke = 0.2, shape = 21) +
    scale_fill_gradient(low = "white", high = "navyblue",
			breaks = c(0, 1), labels = c("Low", "High")) +
    scale_size(range = c(0.1, 3), 
	       limits = c(0, 1),
	       breaks = c(0, .33, .66, 1)) +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 1.1, vjust = 0.5),
	  axis.text.y = element_text(size = 8),
	  panel.grid = element_line(linewidth = .25, color = "grey90"),
	  legend.margin = margin(0, 0, -1, -.25, unit = "lines"),
	  legend.title = element_text(size = 9),
	  legend.text = element_text(size = 8, margin = margin(l = 0, r= 0, unit = "lines")),
	  plot.margin = margin(0, 0, 0, 0)
	  ) +
    guides(fill = guide_colorbar(order = 1, position = "right",
				 barwidth = .25, barheight = 5),
	   size = guide_legend(position = "top")) +
    labs(x = NULL, y = NULL, fill = "Express:", size = "Proportion:")


fig_d_title <- 
    ggdraw() + 
    draw_label(
	       "Marker genes' RNA expression",
	       x = 0,
	       size = 9,
	       hjust = 0
	       ) +
    theme(text = element_text(size = 9),
	  plot.margin = margin(l = 1.25, unit = "lines"))

fig_d <- plot_grid(fig_d_title, dotplot_markers, ncol = 1, rel_heights = c(.1, 1))



# Figure E ####################################################################
gating_subtypes <-
    bcells@meta.data |>
    as_tibble(rownames = "barcode") |>
    select(barcode, starts_with("is.pure_")) |>
    pivot_longer(-barcode, names_to = "subtype") |>
    mutate(subtype = str_remove(subtype, "is.pure_"),
	   value = recode(value, "Impure" = "No", "Pure" = "Yes"))

gating_scores <- 
    bcells@meta.data |>
    as_tibble(rownames = "barcode") |>
    select(barcode, ends_with("UCell")) |>
    pivot_longer(-barcode, names_to = "subtype", values_to = "score") |>
    mutate(subtype = str_remove(subtype, "_UCell"),
	   subtype = recode(subtype, "memory" = "mem", "prolif." = "prolif"))

gating_df <-
    left_join(gating_subtypes, gating_scores, join_by(barcode, subtype)) |>
    group_by(subtype) |>
    mutate(value = case_when(subtype %in% c("naive", "plasma") & value == "Yes" & score >= quantile(score, .8) ~ "Yes",
			     subtype %in% c("naive", "plasma") & value == "Yes" & score < quantile(score, .8) ~ "No",
			     .default = value)) |>
    ungroup() |>
    left_join(select(umap_df, barcode, umap_1, umap_2), join_by(barcode))

subsets_df <-
    gating_df |>
    left_join(select(umap_df, barcode, hto), join_by(barcode)) |>
    select(barcode, subtype, hto, value) |>
    group_by(hto, subtype) |>
    summarise(pct = mean(value == "Yes")) |>
    ungroup() |>
    mutate(subtype = recode(subtype, 
			    "abc" = "CD11c<sup>+</sup>",
			    "mem" = "IgD<sup>–</sup>  CD27<sup>+</sup>",
			    "naive" = "IgD<sup>+</sup>  CD27<sup>–</sup>",
			    "plasma" = "*XBP1*<sup>high</sup>",
			    "prolif" = "*MKI67*<sup>+</sup>")) |>
    arrange(hto, desc(pct)) |>
    mutate(subtype = fct_inorder(subtype))

proportions_plot <- 
    ggplot(subsets_df, aes(x = hto, y = pct)) +
    geom_col(aes(fill = hto)) +
    scale_fill_manual(values = stim_colors) +
    facet_wrap(~subtype, nrow = 1, scales = "free_y") +
    ggh4x::facetted_pos_scales(y = list(subtype == "IgD<sup>+</sup>  CD27<sup>–</sup>" ~ scale_y_continuous(breaks = c(0, 0.4, 0.8)),
					subtype == "IgD<sup>–</sup>  CD27<sup>+</sup>" ~ scale_y_continuous(breaks = c(0, 0.1, 0.2)),
					subtype == "CD11c<sup>+</sup>" ~ scale_y_continuous(limits = c(0, .2), breaks = c(0, 0.1, 0.2)),
					subtype == "*XBP1*<sup>high</sup>" ~ scale_y_continuous(breaks = c(0, 0.25, 0.5)),
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


fig_e_title <- 
    ggdraw() + 
    draw_label(
	       "Proportion of cell subsets in each condition",
	       x = 0,
	       size = 9,
	       hjust = 0
	       ) +
    theme(text = element_text(size = 9),
	  plot.margin = margin(l = 1.25, unit = "lines"))

fig_e <- plot_grid(fig_e_title, proportions_plot, ncol = 1, rel_heights = c(.1, 1),
		   labels = "e", label_size = 12)



# Figure F ####################################################################
disease_genes <- 
    features_df |>
    filter(gene_name %in% c("IRF8", "STAT1", "BANK1", "BLK", "IL4R", "BATF", "MYC"), 
	   phenotype == "Gene Expression") |>
    select(gene_name, gene_id)

dis_genes_expression <-
    LayerData(bcells, assay = "RNA", layer = "data")[disease_genes$gene_id, ] |>
    as_tibble(rownames = "gene_id") |>
    left_join(disease_genes, join_by(gene_id)) |>
    pivot_longer(-c(gene_id, gene_name), names_to = "barcode") |>
    left_join(select(umap_df, barcode, umap_1, umap_2), join_by(barcode)) |>
    group_by(gene_id, gene_name) |>
    mutate(value = clip_expression(value)) |>
    ungroup()

umaps_disease <-
    dis_genes_expression |>
    group_split(gene_name) |>
    map(~ggplot(data = ., aes(x = umap_1, y = umap_2)) +
	geom_point(aes(color = value), size = .2, stroke = 0) +
	scale_color_gradient(low = "Snow", high = "black") +
	facet_wrap(~gene_name) +
	theme_minimal() +
	theme(axis.title = element_blank(),
	      axis.text = element_blank(),
	      legend.text = element_blank(),
	      legend.margin = margin(0, 0, 0, -0.5, "lines"),
	      panel.grid = element_blank(),
	      strip.text = element_text(size = 8, margin = margin(0, 0, 0, 0)),
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
	       "RNA expression of selected disease risk genes",
	       x = 0,
	       size = 9,
	       hjust = 0
	       ) +
    theme(text = element_text(size = 9),
	  plot.margin = margin(l = 1.25, unit = "lines"))

fig_f <- plot_grid(fig_f_title, umaps_disease, ncol = 1, rel_heights = c(.1, 1),
		   labels = "f", label_size = 12)





# Final Figure ################################################################
ggsave("./fig3.png",
       plot_grid(
		 plot_grid(fig_a, NULL, fig_b, nrow = 1, rel_widths = c(.95, .05, 1), 
			   labels = c("a", "", "b"), label_size = 12),
		 NULL, 
		 plot_grid(fig_c, fig_d, nrow = 1, rel_widths = c(.5, 1), 
			   labels = c("c", "d"), label_size = 12, align = "none"),
		 NULL,
		 fig_e,
		 NULL,
		 fig_f,
		 ncol = 1, rel_heights = c(1, 0.05, 1.1, 0.025, 0.4, 0.05, 0.5)) +
       theme(plot.background = element_rect(color = "white", fill = "white")),
       width = 6.5, height = 8, dpi = 600)

