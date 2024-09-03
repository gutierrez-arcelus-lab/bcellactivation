# Set large RAM limit for R
unix::rlimit_as(1e12)

# single-cell data analysis
library(Seurat)
library(SeuratDisk)
library(harmony)

# Data wrangling
library(tidyverse)

# Plotting
library(scales)
library(tidytext)
library(cowplot)
library(patchwork)

# Colors
stim_colors <- 
    "../figure_colors.txt" |>
    read_tsv(col_names = c("stim", "color")) |>
    mutate(stim = sub("_", " ", stim),
	   stim = paste0(stim, "h")) |>
    deframe()

cellranger_dir_1984 <- 
    file.path("/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/B_cells_citeseq",
	      "SN0268787/broad/hptmp/curtism/bwh10x/KW10456_Maria",
	      "221101_10X_KW10456_bcl/cellranger-6.1.1/GRCh38/BRI-1984/outs",
	      "filtered_feature_bc_matrix")

features_df <- file.path(cellranger_dir_1984, "features.tsv.gz") |>
    read_tsv(col_names = c("gene_id", "gene_name", "phenotype"))

genes_df <- features_df |>
    filter(phenotype == "Gene Expression") |>
    select(gene_id, gene_name)

mt_ribo_genes <- 
    genes_df |>
    filter(grepl("^MT-|^MRPS|^MRPL|^RPS|^RPL", gene_name))

# Import Seurat object
bcells <- read_rds("./data/seurat_qced.rds")

# Scale and run PCA
bcells <- bcells |>
    FindVariableFeatures() |>
    ScaleData(vars.to.regress = c("nCount_RNA", "percent_mt")) |>
    RunPCA()

# Run Harmony correcting for batch
set.seed(1L)
bcells <- bcells |>
    RunHarmony(group.by.vars = c("orig.ident", "donor_id"),
	       max.iter.harmony = 30,
	       reduction.save = "harmony") |>
    FindNeighbors(dims = 1:30, reduction = "harmony") |>
    FindClusters(resolution = 0.5) |>
    RunUMAP(reduction = "harmony", 
	    dims = 1:30,
	    seed.use = 1L,
	    reduction.name = "umap")

umap_df <- 
    Embeddings(bcells, "umap") |>
    as_tibble(rownames = "barcode") |>
    left_join(as_tibble(bcells@meta.data, rownames = "barcode"), join_by(barcode)) |>
    select(barcode, lib = orig.ident, donor_id, hto = dmm_hto_call, cluster = RNA_snn_res.0.5,
	   UMAP_1, UMAP_2) |>
    sample_frac(1L) |>
    mutate(barcode = fct_inorder(barcode),
	   hto = factor(hto, levels = names(stim_colors)),
	   cluster = factor(cluster, levels = sort(as.integer(levels(cluster)))))

write_tsv(umap_df, "./data/umap_df.tsv")

# Marker genes
Idents(bcells) <- "RNA_snn_res.0.5"

cluster_markers <- 
    FindAllMarkers(bcells, 
                   only.pos = TRUE,
                   min.pct = 0.1,
                   logfc.threshold = .5) |>
    as_tibble() |>
    left_join(genes_df, join_by(gene == gene_id))

write_tsv(cluster_markers, "./data/cluster_markers.tsv")



# Plots
umap_batch <-
    ggplot(umap_df, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(color = lib), size = .1) +
    scale_color_manual(values = c("#C5000B", "midnightblue", "#FFB022")) +
    theme_minimal() +
    theme(axis.title = element_blank(),
	  axis.text = element_blank(),
	  strip.text.y = element_blank(),
	  axis.ticks = element_blank(),
	  panel.grid = element_blank(),
	  plot.background = element_rect(fill = "white", color = "white")) +
    guides(color = guide_legend(title = "Batch:", 
				override.aes = list(size = 4, alpha = 1)))

umap_donor <-
    ggplot(umap_df, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(color = donor_id), size = .1) +
    scale_color_manual(values = RColorBrewer::brewer.pal(n = 5, "Set1")) +
    theme_minimal() +
    theme(axis.title = element_blank(),
	  axis.text = element_blank(),
	  strip.text.y = element_blank(),
	  axis.ticks = element_blank(),
	  panel.grid = element_blank(),
	  plot.background = element_rect(fill = "white", color = "white")) +
    guides(color = guide_legend(title = "Donor ID:", 
				override.aes = list(size = 4, alpha = 1)))

umap_stim <-
    ggplot(umap_df, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(color = hto), size = .1) +
    scale_color_manual(values = stim_colors) +
    theme_minimal() +
    theme(axis.title = element_blank(),
	  axis.text = element_blank(),
	  strip.text.y = element_blank(),
	  axis.ticks = element_blank(),
	  panel.grid = element_blank(),
	  plot.background = element_rect(fill = "white", color = "white")) +
    guides(color = guide_legend(title = "Stim:", 
				override.aes = list(size = 4, alpha = 1)))

# Clusters
umap_clust <- 
    ggplot(umap_df, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(color = cluster), size = .1) +
    scale_color_manual(values = pals::kelly(n = length(unique(umap_df$cluster)) + 1 )[-1] ) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
	  plot.title = element_text(hjust = .5),
	  plot.background = element_rect(fill = "white", color = "white")
	  ) +
    guides(color = guide_legend(title = "Cluster:", ncol = 2, 
				override.aes = list(size = 4, alpha = 1)))

mki_gene_df <- features_df |>
    filter(gene_name == "MKI67") |>
    select(gene_id, gene_name)

mki_gene_quant <- bcells@assays$RNA@data |> 
    {function(x) x[mki_gene_df$gene_id, , drop = FALSE]}() |>
    as_tibble(rownames = "gene_id") |>
    left_join(mki_gene_df, by = "gene_id") |>
    pivot_longer(-c("gene_id", "gene_name"), 
                 names_to = "barcode", values_to = "gene_exp")

umap_ki67 <- umap_df |> 
    left_join(mki_gene_quant, by = "barcode") |>
    select(barcode, gene_exp, UMAP_1, UMAP_2) |>
    ggplot(aes(UMAP_1, UMAP_2)) +
    geom_point(aes(color = gene_exp), size = .1) +
    scale_color_gradient(low = "grey80", high = "midnightblue", 
			 guide = guide_colorbar(direction = "vertical",
						barwidth = .5)) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
	  plot.background = element_rect(fill = "white", color = "white")
	  ) +
    labs(color = "MKI67:")

umap_df_features <- bcells@meta.data |> 
    as_tibble(rownames = "barcode") |>
    left_join(umap_df, join_by(barcode)) |>
    select(barcode, n_genes = nFeature_RNA, percent_mt, UMAP_1, UMAP_2)

umap_genes <- 
    ggplot(umap_df_features, 
	   aes(UMAP_1, UMAP_2)) +
    geom_point(aes(color = n_genes), size = .1) +
    scale_color_gradient(low = "grey80", high = "midnightblue", 
			 guide = guide_colorbar(direction = "vertical",
						barwidth = .5)) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
	  plot.background = element_rect(fill = "white", color = "white")
	  ) +
    labs(color = "Genes:")

umap_mito <- 
    ggplot(umap_df_features, 
	   aes(UMAP_1, UMAP_2)) +
    geom_point(aes(color = percent_mt), size = .1) +
    scale_color_gradient(low = "grey80", high = "midnightblue", 
			 guide = guide_colorbar(direction = "vertical",
						barwidth = .5)) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
	  plot.background = element_rect(fill = "white", color = "white")
	  ) +
    labs(color = "% Mito:")

ggsave("./plots/umap.png", 
       umap_batch + umap_stim + umap_clust + umap_ki67 + umap_genes + umap_mito + plot_layout(ncol = 2),
       width = 8, height = 7, dpi = 600) 



# Marker genes
top_markers <- cluster_markers |>
    as_tibble() |>
    group_by(cluster) |>
    top_n(10, avg_log2FC) |>
    ungroup() |>
    select(cluster, gene_name, gene_id = gene, avg_log2FC)
  
cluster_cell <- Idents(bcells) |>
    enframe(name = "barcode", value = "cluster_cell") |>
    mutate(cluster_cell = factor(cluster_cell, levels = levels(top_markers$cluster)))
  
cell_gene_expr <- bcells@assays$RNA@data |>
    {function(x) x[unique(top_markers$gene_id), ]}() |>
    as_tibble(rownames = "gene_id") |>
    pivot_longer(-gene_id, names_to = "barcode", values_to = "logexpr")
  
top_marker_expr <- cell_gene_expr |>
    left_join(cluster_cell, join_by(barcode)) |>
    inner_join(top_markers, join_by(gene_id), relationship = "many-to-many")
  
top_marker_perc_exp <- top_marker_expr |>
    group_by(cluster = cluster_cell, gene_name) |>
    summarise(prop_expr = mean(logexpr > 0)) |>
    ungroup()
  
top_marker_avg_exp <- top_marker_expr |>
    filter(logexpr > 0) |>
    group_by(cluster = cluster_cell, gene_name) |>
    summarise(scaled_expr = mean(logexpr)) |>
    ungroup()
  
top_marker_summary <- 
    left_join(top_marker_perc_exp, top_marker_avg_exp, join_by(cluster, gene_name)) |>
    left_join(select(top_markers, cluster_top = cluster, gene_name, avg_log2FC), 
	      join_by(gene_name),
	      relationship = "many-to-many") |>
    arrange(cluster_top, avg_log2FC) |>
    mutate(gene_name = fct_inorder(gene_name)) 
  
cluster_markers_plot <-
    ggplot(top_marker_summary, aes(cluster, gene_name)) +
        geom_point(aes(size = prop_expr, fill = scaled_expr), 
                   color = "black", shape = 21) +
        scale_size(range = c(0.1, 4), labels = scales::percent) +
        scale_fill_viridis_c(option = "magma") +
        facet_wrap(~cluster_top, scales = "free", ncol = 3) +
        theme_bw() +
        theme(axis.line = element_blank(),
              axis.text.x = element_text(size = 8, angle = 45, hjust = 1, vjust = 1),
              axis.text.y = element_text(size = 8),
              panel.grid = element_line(color = "grey96"),
              panel.border = element_blank(),
              legend.position = "top") +
        labs(x = NULL, y = NULL, 
             fill = "Scaled\nExpression", 
             size = "% of cells") +
        guides(fill = guide_colorbar(barheight = .5))

ggsave("./plots/cluster_markers.png", 
       cluster_markers_plot,
       height = 10, width = 8, dpi = 600)

# ADT
adt_markers <- 
    FindAllMarkers(bcells, 
		   assay = "ADT",
                   only.pos = TRUE,
                   min.pct = 0.1,
                   logfc.threshold = .5) |>
    as_tibble() 

clip_expression <- function(values) { 
    
    q_clip <- quantile(values, 0.99)

    case_when(values >= q_clip ~ q_clip,
	      TRUE ~ values)
}

adt_df <- 
    bcells@assays$ADT@data |> 
    t() |>
    as_tibble(rownames = "barcode") |>
    pivot_longer(-barcode, names_to = "adt") |>
    left_join(umap_df, join_by(barcode)) |>
    select(barcode, UMAP_1, UMAP_2, adt, value) |>
    group_by(adt) |>
    mutate(value = clip_expression(value)) |>
    ungroup() |>
    sample_frac(1L) |>
    mutate(adt = factor(adt, levels = rownames(bcells@assays$ADT@data)),
	   barcode = fct_inorder(barcode))
    
adt_plot_list <- 
    adt_df |>
    group_split(adt) |>
    map(~ggplot(data = ., aes(x = UMAP_1, y = UMAP_2)) +
	geom_point(aes(color = value), size = .2, stroke = 0) +
	scale_color_gradient(low = "grey80", high = "darkred") +
	facet_wrap(~adt, ncol = 1) +
	theme_minimal() +
	theme(axis.title = element_blank(),
	      axis.text = element_blank(),
	      legend.text = element_text(size = 7, margin = margin(l = 0.1, unit = "lines")),
	      legend.margin = margin(0, 0, 0, -0.5, "lines"),
	      panel.grid = element_blank(),
	      strip.text = element_text(size = 6, margin = margin(0, 0, 0, 0)),
	      strip.clip = "off",
	      panel.spacing = unit(0, "lines"),
	      plot.margin = margin(0, 0, 0, 0)
	      ) +
	labs(color = NULL) +
	guides(color = guide_colorbar(barwidth = .2, barheight = 2))
    )


ggsave("./plots/adt_3.png",
       plot_grid(plotlist = adt_plot_list[97:137], ncol = 6) +
	   theme(plot.margin = margin(1, .5, 1, .5, unit = "in"),
		 plot.background = element_rect(color = "white", fill = "white")),
       width = 8.5, height = 11, dpi = 600)


###
FindMarkers(bcells, 
	    ident.1 = "10",
	    ident.2 = "9",
	    only.pos = TRUE,
	    min.pct = 0.1,
	    logfc.threshold = .5) |>
as_tibble(rownames = "gene_id") |>
left_join(genes_df, join_by(gene_id))


bcells@reductions$harmony@feature.loadings |>
    as_tibble(rownames = "gene_id") |>
    left_join(genes_df) |>
    pivot_longer(-c(gene_id, gene_name), names_to = "dim") |>
    group_by(dim) |>
    top_n(4, value) |>
    ungroup() |>
    mutate(dim = str_remove(dim, "harmony_"),
	   dim = as.integer(dim)) |>
    arrange(dim, desc(value)) |>
    slice(41:60)



###




