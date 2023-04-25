# Set large RAM limit for R
unix::rlimit_as(1e12)

# single-cell data analysis
library(Seurat)
library(harmony)
devtools::load_all("/lab-share/IM-Gutierrez-e2/Public/vitor/sclibr")

# Data wrangling
library(dplyr)
library(forcats)
library(tidyr)
library(purrr)
library(readr)
library(tibble)
library(ggplot2)

# Plotting
library(scales)
library(tidytext)
library(ggridges)
library(RColorBrewer)
library(scico)
library(ggsci)
library(ggthemes)
library(ggnewscale)
library(cowplot)


# Colors
stim_order <- 
    c("Day 0", 
      sprintf("IL4 %sh", c(24, 72)),
      sprintf("BCR %sh", c(24, 72)),
      sprintf("TLR7 %sh", c(24, 72)),
      sprintf("DN2 %sh", c(24, 72))
    )
      
stim_colors <- 
    c("grey80",
      "grey50", "black",
      brewer.pal(n = 9, "Blues")[c(3, 8)],
      brewer.pal(n = 9, "Greens")[c(3, 8)],
      paste0("tomato", c(2, 4))
      )

names(stim_colors) <- stim_order


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

# Import Seurat objects
lib1984 <- read_rds("./data/seurat_1984_qced.rds")
lib1988 <- read_rds("./data/seurat_1988_qced.rds")
lib1990 <- read_rds("./data/seurat_1990_qced.rds")

DefaultAssay(lib1984) <- DefaultAssay(lib1988) <- DefaultAssay(lib1990) <- "RNA"

# Merge datasets and run PCA
bcells <- 
    merge(lib1984, y = c(lib1988, lib1990), add.cell.ids = c("1984", "1988", "1990")) |>
    {function(x) ScaleData(x, features = rownames(x))}() |>
    FindVariableFeatures() |>
    {function(x) RunPCA(x, features = VariableFeatures(x))}()

# Visualize PCA
variable_genes_merge <- VariableFeatures(bcells)

total_variance <- bcells@assays$RNA@scale.data[variable_genes_merge, ] |>
    matrixStats::rowVars() |>
    sum()

sdev_plot <- tibble(sdev = bcells@reductions$pca@stdev) |>
    rownames_to_column("pc") |>
    mutate(pc = fct_inorder(pc)) |>
    ggplot(aes(pc, sdev)) +
    geom_line(aes(group = 1)) +
    scale_x_discrete(breaks = c(1, seq(5, 50, 5))) +
    scale_y_continuous(sec.axis = sec_axis(trans = ~ (.^2) / total_variance,
					   labels = percent,
					   name = "% variance explained")) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    labs(x = "Principal component", y = "Standard deviation")

ggsave("./plots/pca_sdev.png", sdev_plot, width = 3.5, height = 2.5)

# Run Harmony
# Correct for batch
bcells <- bcells |>
    RunHarmony(group.by.vars = "orig.ident", 
	       max.iter.harmony = 30,
	       reduction.save = "harmony.1")

# Correct for batch and condition
bcells <- bcells |>
    RunHarmony(group.by.vars = c("orig.ident", "hto"), 
	       max.iter.harmony = 30,
	       reduction.save = "harmony.2")

# UMAP
bcells <- bcells |>
    RunUMAP(reduction = "pca", 
	    dims = 1:35, 
	    seed.use = 1,
	    reduction.name = "umap.merged")

bcells <- bcells |>
    RunUMAP(reduction = "harmony.1", 
	    dims = 1:35,
	    seed.use = 1,
	    reduction.name = "umap.harmony.1")

bcells <- bcells |>
    RunUMAP(reduction = "harmony.2", 
	    dims = 1:35,
	    seed.use = 1,
	    reduction.name = "umap.harmony.2")

umap_df <- 
    list("Merge" = Embeddings(bcells, "umap.merged"),
	 "Harmony (batch)" = Embeddings(bcells, "umap.harmony.1"),
	 "Harmony (batch + stim)" = Embeddings(bcells, "umap.harmony.2")) |>
    map_dfr(~as_tibble(., rownames = "barcode") |> 
		select(barcode, UMAP_1 = 2, UMAP_2 = 3),
	    .id = "reduction") |>
    left_join(as_tibble(bcells@meta.data, rownames = "barcode")) |>
    select(reduction, barcode, UMAP_1, UMAP_2, dataset = orig.ident, hto) |>
    mutate(reduction = fct_inorder(reduction),
	   dataset = paste0("BRI-", dataset)) |>
    sample_frac(1L) |>
    mutate(barcode = fct_inorder(barcode))



umap <- umap_df |>
    pivot_longer(dataset:hto, names_to = "variable") |>
    ggplot(aes(UMAP_1, UMAP_2)) +
    geom_point(data = ~subset(., variable == "dataset"),
	       aes(color = value), 
	       size = .1, alpha = .25) +
    scale_color_manual(values = c("black", "turquoise2", "deeppink3")) +
    guides(color = guide_legend(title = "dataset", 
				order = 1, 
				override.aes = list(size = 2, alpha = 1))) +
    new_scale_color() +
    geom_point(data = ~subset(., variable == "hto") |>
		mutate(value = factor(value, levels = stim_order)), 
	       aes(color = value), 
	       size = .1, alpha = .25) +
    scale_color_manual(values = stim_colors) +
    guides(color = guide_legend(title = "stim",
				order = 2, 
				override.aes = list(size = 2, alpha = 1))) +
    facet_grid(variable~reduction) +
    theme_bw() +
    theme(axis.title = element_blank(),
	  axis.text = element_blank(),
	  strip.text.y = element_blank(),
	  axis.ticks = element_blank(),
	  panel.grid = element_blank(),
	  plot.background = element_rect(fill = "white", color = "white"))

ggsave("./plots/umap.png", umap, width = 8, height = 5, dpi = 600) 



# Clusters
bcells <- bcells |>
    FindNeighbors(dims = 1:35, reduction = "harmony.1", nn.eps = .5) |>
    FindClusters(resolution = 0.4)

cluster_df <- bcells@meta.data |>
    as_tibble(rownames = "barcode") |>
    select(barcode, cluster = RNA_snn_res.0.4)

umap_df_clu <- umap_df |>
    filter(reduction == "Harmony (batch)") |>
    left_join(cluster_df, join_by(barcode)) |>
    select(barcode, cluster, UMAP_1, UMAP_2)


# Plot clusters and features
sc1 <- range(umap_df_clu$UMAP_1) |>
    {function(x) length(x[1]:x[2])}() |>
    {function(x) x/6}()

sc2 <- range(umap_df_clu$UMAP_2) |>
    {function(x) length(x[1]:x[2])}() |>
    {function(x) x/6}()

x1 <- min(umap_df_clu$UMAP_1)
x2 <- x1 + sc1
y1 <- min(umap_df_clu$UMAP_2)
y2 <- y1 + sc2 

x1 <- x1 + 1
x2 <- x2 + 1
y1 <- y1 + 1
y2 <- y2 + 1

cluster_labs <- umap_df_clu |> 
    group_by(cluster) |>
    summarise_at(vars(UMAP_1, UMAP_2), median) |>
    ungroup()

clust_colors <- 
    c("black", "grey", "purple2", "goldenrod2", 
      "magenta3", "yellow", "green4", "royalblue",
      pal_npg()(9)) 

umap_clust <- ggplot(umap_df_clu, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(color = cluster), size = .25) +
    geom_label(data = cluster_labs, aes(label = cluster), alpha = .5, size = 3) +
    scale_color_manual(values = clust_colors) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
	  plot.title = element_text(hjust = .5),
	  legend.position = "none") +
    labs(title = "Clusters")

mki_gene_df <- features_df |>
    filter(gene_name == "MKI67") |>
    select(gene_id, gene_name)

mki_gene_quant <- bcells@assays$RNA@data |> 
    {function(x) x[mki_gene_df$gene_id, ,drop = FALSE]}() |>
    as_tibble(rownames = "gene_id") |>
    left_join(mki_gene_df, by = "gene_id") |>
    pivot_longer(-c("gene_id", "gene_name"), 
                 names_to = "barcode", values_to = "gene_exp")

umap_ki67 <- umap_df_clu |> 
    left_join(mki_gene_quant, by = "barcode") |>
    select(barcode, gene_exp, UMAP_1, UMAP_2) |>
    ggplot(aes(UMAP_1, UMAP_2)) +
    geom_point(aes(color = gene_exp, fill = gene_exp), size = .1, shape = 19) +
    scale_color_gradient(low = "beige", high = "tomato4", 
			 guide = "none") +
    scale_fill_gradient(low = "beige", high = "tomato4", 
			guide = guide_colorbar(barwidth = .5)) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
	  legend.position = c(.9, .5),
	  plot.title = element_text(hjust = .5)) +
    labs(fill = NULL, title = "KI67")


umap_df_features <- bcells@meta.data |> 
    as_tibble(rownames = "barcode") |>
    left_join(umap_df_clu, join_by(barcode)) |>
    select(barcode, n_genes = nFeature_RNA, percent_mt, UMAP_1, UMAP_2)

umap_genes <- 
    ggplot(umap_df_features, 
	   aes(UMAP_1, UMAP_2)) +
    geom_point(aes(color = n_genes, fill = n_genes), size = .1, shape = 19) +
    scale_color_gradient(low = "beige", high = "tomato4", 
			 guide = "none") +
    scale_fill_gradient(low = "beige", high = "tomato4", 
			guide = guide_colorbar(barwidth = .5)) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
	  legend.position = c(.9, .5),
	  plot.title = element_text(hjust = .5)) +
    labs(fill = NULL, title = "Number of genes")

umap_mito <- 
    ggplot(umap_df_features, 
	   aes(UMAP_1, UMAP_2)) +
    geom_point(aes(color = percent_mt, fill = percent_mt), size = .1, shape = 19) +
    scale_color_gradient(low = "beige", high = "tomato4", 
			 guide = "none") +
    scale_fill_gradient(low = "beige", high = "tomato4", 
			guide = guide_colorbar(barwidth = .5)) +
    geom_segment(x = x1, xend = x2, y = y1, yend = y1,
		 arrow = arrow(length = unit(5, "pt")),
		 show.legend = FALSE) +
    geom_segment(x = x1, xend = x1, y = y1, yend = y2,
		 arrow = arrow(length = unit(5, "pt")), 
		 show.legend = FALSE) +
    annotate("text", label = "UMAP 1", size = rel(3),
	     x = x1, y = y1, hjust = 0, vjust = 1.5) + 
    annotate("text", label = "UMAP 2", size = rel(3),
	     x = x1, y = y1, angle = 90, hjust = 0, vjust = -.5) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
	  legend.position = c(.9, .5),
          axis.text = element_blank(),
	  plot.title = element_text(hjust = .5)) +
    labs(fill = NULL, title = "% Mito")
    

umaps <- 
    plot_grid(NULL,
	      umap_clust, 
	      umap_ki67,
	      umap_genes, 
	      umap_mito,
	      ncol = 1, 
	      rel_heights = c(.2, 1, 1, 1, 1)) +
    theme(plot.background = element_rect(fill = "white", color = "white"))


# Marker genes
Idents(bcells) <- "RNA_snn_res.0.4"

cluster_markers <- 
    FindAllMarkers(bcells, 
                   only.pos = TRUE,
                   min.pct = 0.1,
                   logfc.threshold = .5) |>
    as_tibble() |>
    left_join(filter(features_df, phenotype == "Gene Expression"), join_by(gene == gene_id)) |>
    select(-phenotype)

cluster_markers_plot <- plot_markers(bcells, cluster_markers)

ggsave("./plots/cluster_markers.png", 
       plot_grid(umaps, cluster_markers_plot,
		 nrow = 1, rel_widths = c(.33, 1),
		 labels = c("A)", "B)")),
       height = 10, width = 10, dpi = 600)


# B cell subset marker genes
#clip_expression <- function(values) { 
#    
#    q01 <- quantile(values, 0.01)
#    q99 <- quantile(values, 0.99)
#
#    case_when(values <= q01 ~ q01,
#	      values >= q99 ~ q99,
#	      TRUE ~ values)
#}

umap_df2 <- umap_df |>
    filter(reduction == "Harmony (batch + stim)") |>
    select(barcode, hto, UMAP_1, UMAP_2) |>
    mutate(hto = factor(hto, levels = stim_order))

bcell_genes <- 
    c("CD27" = "CD27",
      "CD69" = "CD69",
      "CD83" = "CD83",
      "CD23" = "FCER2",
      "IL4R" = "IL4R",
      "CD11c" = "ITGAX",
      "IgJ" = "JCHAIN",
      "MZB1" = "MZB1",
      "T-bet" = "TBX21",
      "TCL1A" = "TCL1A",
      "TACI" = "TNFRSF13B",
      "BCMA (CD269)" = "TNFRSF17",
      "CD38" = "CD38",
      "CD138" = "SDC1",
      "BANK1" = "BANK1")

bcell_genes_df <- features_df |>
    filter(gene_name %in% bcell_genes, phenotype == "Gene Expression") |>
    select(gene_id, gene_name)

bcell_gene_quant <- bcells@assays$RNA@data |> 
    {function(x) x[bcell_genes_df$gene_id, , drop = FALSE]}() |>
    as_tibble(rownames = "gene_id") |>
    left_join(bcell_genes_df, by = "gene_id") |>
    pivot_longer(-c("gene_id", "gene_name"), 
                 names_to = "barcode", values_to = "gene_exp") |>
    group_by(gene_id) |>
    mutate(gene_exp = clip_expression(gene_exp)) |>
    ungroup() |>
    select(barcode, gene_id, gene_name, gene_exp)

bcell_genes_plot_list <- bcell_gene_quant |>
    left_join(umap_df2) |>
    {function(x) split(x, x$gene_name)}() |>
    map(~ggplot(data = ., aes(UMAP_1, UMAP_2, color = gene_exp)) +
	    geom_point(size = .1) +
	    scale_color_gradient(low = "lightyellow2", high = "tomato4",
				 guide = guide_colorbar(barwidth = .5)) +
	    facet_wrap(~gene_name) +
	    theme_bw() +
	    theme(panel.grid = element_blank(),
		  panel.border = element_blank(),
		  axis.title = element_blank(),
		  axis.text = element_blank(),
		  strip.background = element_blank(),
		  axis.line = element_blank(),
		  axis.ticks = element_blank(),
		  strip.text = element_text(size = 10, face = "bold"),
		  legend.position = "none"))

umap_stim_int <- umap_df2 |>
    mutate(lab = "Stim") |>
    ggplot(aes(UMAP_1, UMAP_2)) +
    geom_point(aes(color = hto), size = .1, alpha = 1) +
    scale_color_manual(values = stim_colors) +
    facet_wrap(~lab) +
    theme_bw() +
	    theme_bw() +
	    theme(panel.grid = element_blank(),
		  panel.border = element_blank(),
		  axis.title = element_blank(),
		  axis.text = element_blank(),
		  strip.background = element_blank(),
		  axis.line = element_blank(),
		  axis.ticks = element_blank(),
		  strip.text = element_text(size = 10, face = "bold"),
		  legend.position = "none")

bcell_genes_plot <- plot_grid(plotlist = c(list(umap_stim_int), bcell_genes_plot_list)) +
    theme(plot.background = element_rect(fill = "white", color = "white"))

ggsave("./plots/bcell_genes.png", bcell_genes_plot, width = 10, height = 10, dpi = 600)


# ADT
bcell_prots <- 
    c("IGHD" = "IgD", 
      "IGHM" = "IgM", 
      "CR2" = "CD21", 
      "FCER2" = "CD23",
      "CD27" = "CD27",
      "CD38" = "CD38",
      "CD45RA" = "CD45RA",
      "CD69" = "CD69", 
      "CD86" = "CD86",
      "CD99" = "CD99",
      "ITGAL" = "CD11a",
      "ITGAM" = "CD11b",
      "ITGAX" = "CD11c",
      "CXCR5" = "CD185-or-CXCR5",
      "PDCD1" = "CD279",
      "FCGR2A" = "CD32",
      "CTLA4" = "CD152-or-CTLA-4",
      "TNFRSF13B" = "CD267-or-TACI",
      "NT5E" = "CD73-or-Ecto-5-nucleotidase",
      "GGT1" = "CD224")

bcell_prot_quant <- bcells@assays$ADT@data |> 
    {function(x) x[bcell_prots, , drop = FALSE]}() |>
    as_tibble(rownames = "gene_id") |>
    pivot_longer(-gene_id, names_to = "barcode", values_to = "prot_exp") |>
    group_by(gene_id) |>
    mutate(prot_exp = clip_expression(prot_exp)) |>
    ungroup() |>
    select(barcode, gene_id, prot_exp)

bcell_prot_plot_list <- bcell_prot_quant |>
    left_join(umap_df2) |>
    {function(x) split(x, x$gene_id)}() |>
    map(~ggplot(data = ., aes(UMAP_1, UMAP_2, color = prot_exp)) +
	    geom_point(size = .1) +
	    scale_color_scico(palette = "lajolla",
			      labels = function(x) str_pad(x, 3),
			      guide = guide_colorbar(barwidth = .5)) +
	    facet_wrap(~gene_id) +
	    theme_bw() +
	    theme(panel.grid = element_blank(),
		  panel.border = element_blank(),
		  axis.title = element_blank(),
		  axis.text = element_blank(),
		  strip.background = element_blank(),
		  axis.line = element_blank(),
		  axis.ticks = element_blank(),
		  strip.text = element_text(size = 10, face = "bold"),
		  legend.position = "none") +
	    labs(color = NULL))

bcell_prot_plot <- plot_grid(plotlist = c(list(umap_stim_int), bcell_prot_plot_list)) +
    theme(plot.background = element_rect(fill = "white", color = "white"))

ggsave("./plots/bcell_prots.png", bcell_prot_plot, width = 10, height = 10, dpi = 600)






## Seurat integration with CCA
#
## Integrate batches
#bcell_objects <- list("1984" = lib1984,
#		      "1988" = lib1988,
#		      "1990" = lib1990)
#
#integration_features <- SelectIntegrationFeatures(object.list = bcell_objects)
#
#anchors <- FindIntegrationAnchors(object.list = bcell_objects,
#				  anchor.features = integration_features,
#				  assay = c("RNA", "RNA", "RNA"),
#				  reduction = "cca",
#				  k.anchor = 20,
#				  dims = 1:35)
#
#bcells_integrated <- 
#    IntegrateData(anchorset = anchors, 
#		  dims = 1:35,
#		  features.to.integrate = rownames(lib1984@assays$RNA@counts))
#
#write_rds(bcells_integrated, "./data/seurat_mgb_integrated.rds")
##bcells_integrated <- read_rds("./data/seurat_mgb_integrated.rds")
#
## PCA on integrated data
#bcells_integrated <- bcells_integrated |>
#    {function(x) ScaleData(x, features = rownames(x))}() |>
#    {function(x) RunPCA(x, features = VariableFeatures(x))}() |>
#    RunUMAP(dims = 1:35, verbose = FALSE)
#
## UMAP
#umap_df_s <- bcells_integrated@reductions$umap@cell.embeddings |>
#    as_tibble(rownames = "barcode") |>
#    left_join(as_tibble(bcells_integrated@meta.data, rownames = "barcode")) |>
#    select(barcode, UMAP_1, UMAP_2, dataset = orig.ident, hto) |>
#    mutate(hto = factor(hto, levels = stim_order))
#
#umap_stims_s <- plot_umap(umap_df_s, hto) +
#    scale_color_manual(values = stim_colors) +
#    theme(legend.position = c(.75, .25),
#          legend.key.height = unit(.75, "lines")) +
#    labs(color = "Stim")
#
#umap_batch_s <- plot_umap(umap_df_s, dataset) + 
#    scale_color_manual(values = c("goldenrod2", "red4", "royalblue")) +
#    theme(legend.position = c(.75, .25)) +
#    labs(color = "batch")
#
#ggsave("./plots/umap_seuratIntegration.png", 
#       plot_grid(umap_batch_s, umap_stims_s, nrow = 1), 
#       width = 8, height = 3)
#
