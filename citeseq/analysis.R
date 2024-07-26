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
library(ggridges)
library(RColorBrewer)
library(scico)
library(ggsci)
library(ggthemes)
library(ggnewscale)
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
    FindNeighbors(dims = 1:30, reduction = "harmony", nn.eps = .5) |>
    FindClusters(resolution = 0.4) |>
    RunUMAP(reduction = "harmony", 
	    dims = 1:30,
	    seed.use = 1L,
	    reduction.name = "umap")

# Plots
umap_df <- 
    Embeddings(bcells, "umap") |>
    as_tibble(rownames = "barcode") |>
    left_join(as_tibble(bcells@meta.data, rownames = "barcode"), join_by(barcode)) |>
    select(barcode, lib = orig.ident, donor_id, hto = dmm_hto_call, UMAP_1, UMAP_2) |>
    sample_frac(1L) |>
    mutate(barcode = fct_inorder(barcode),
	   hto = factor(hto, levels = names(stim_colors)))

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
cluster_df <- bcells@meta.data |>
    as_tibble(rownames = "barcode") |>
    select(barcode, cluster = RNA_snn_res.0.4)

umap_df_clu <- umap_df |>
    left_join(cluster_df, join_by(barcode)) |>
    select(barcode, cluster, UMAP_1, UMAP_2)

umap_clust <- 
    ggplot(umap_df_clu, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(color = cluster), size = .1) +
    scale_color_manual(values = pals::kelly(n = length(unique(cluster_df$cluster)) + 1 )[-1] ) +
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

umap_ki67 <- umap_df_clu |> 
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
    left_join(umap_df_clu, join_by(barcode)) |>
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
plot_markers <- function(seurat_obj, cluster_df) {
  
    top_markers <- cluster_df |>
        as_tibble() |>
        group_by(cluster) |>
        top_n(10, avg_log2FC) |>
        ungroup() |>
        select(cluster, gene_name, gene_id = gene, avg_log2FC)
  
    cluster_cell <- Idents(seurat_obj) |>
        enframe(name = "barcode", value = "cluster_cell") |>
        mutate(cluster_cell = factor(cluster_cell, levels = levels(top_markers$cluster)))
  
    cell_gene_expr <- seurat_obj@assays$RNA@data |>
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
        mutate_at(vars(cluster, cluster_top), factor) |>
        arrange(cluster_top, avg_log2FC) |>
        mutate(gene_name = fct_inorder(gene_name)) 
  
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
}

Idents(bcells) <- "RNA_snn_res.0.4"

cluster_markers <- 
    FindAllMarkers(bcells, 
                   only.pos = TRUE,
                   min.pct = 0.1,
                   logfc.threshold = .5) |>
    as_tibble() |>
    left_join(genes_df, join_by(gene == gene_id))

cluster_markers_plot <- plot_markers(bcells, cluster_markers)

ggsave("./plots/cluster_markers.png", 
       cluster_markers_plot,
       height = 10, width = 8, dpi = 600)






# B cell subset marker genes
clip_expression <- function(values) { 
    
    q01 <- quantile(values, 0.01)
    q99 <- quantile(values, 0.99)

    case_when(values <= q01 ~ q01,
	      values >= q99 ~ q99,
	      TRUE ~ values)
}

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

# Poster
#sle_genes <- 
#    "../sle_variants/paper_data/langefeld_top.tsv" |>
#    read_tsv() |>
#    select(locus) |>
#    separate_rows(locus, sep = "-") |>
#    pull(locus)
#
sle_genes <- 
    "../colocalization/finemap/data/bentham_leadvars.tsv" |>
    read_tsv() |>
    select(locus) |>
    separate_rows(locus, sep = ", ") |>
    pull(locus)

sle_genes <- c("BANK1", "IRF8", "STAT1", "SOCS1")

sle_genes_df <- features_df |>
    filter(gene_name %in% sle_genes, phenotype == "Gene Expression") |>
    select(gene_id, gene_name)

sle_gene_quant <- bcells@assays$RNA@data |> 
    {function(x) x[sle_genes_df$gene_id, , drop = FALSE]}() |>
    as_tibble(rownames = "gene_id") |>
    left_join(sle_genes_df, by = "gene_id") |>
    pivot_longer(-c("gene_id", "gene_name"), 
                 names_to = "barcode", values_to = "gene_exp") |>
    select(barcode, gene_id, gene_name, gene_exp) |>
    sample_frac(1L) |>
    mutate(barcode = fct_inorder(barcode))

sle_genes_plot_list <- 
    sle_gene_quant |>
    left_join(umap_df) |>
    {function(x) split(x, x$gene_name)}() |>
    map(~ggplot(data = ., aes(UMAP_1, UMAP_2, color = gene_exp)) +
	    geom_point(size = .1) +
	    scale_color_gradient(low = "grey90", high = "midnightblue",
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
		  strip.text = element_text(size = 10, face = "bold")))

sle_genes_plot <- plot_grid(plotlist = sle_genes_plot_list) +
    theme(plot.background = element_rect(fill = "white", color = "white"))

ggsave("./plots/sle_genes.png", sle_genes_plot, width = 9, height = 8, dpi = 600)





umap_ki67 <- mki_gene_quant |>
    left_join(umap_df3) |>
    select(barcode, gene_id, gene_name, gene_exp, UMAP_1, UMAP_2) |>
    sample_frac(1) |>
    mutate(barcode = fct_inorder(barcode)) |>
    ggplot(aes(UMAP_1, UMAP_2)) +
    geom_point(aes(color = gene_exp, fill = gene_exp), size = .2, shape = 19) +
    scale_color_gradient(low = "grey90", high = "midnightblue", 
			 guide = "none") +
    scale_fill_gradient(low = "grey90", high = "midnightblue", 
			guide = guide_colorbar(barwidth = .5)) +
    facet_wrap(~gene_name) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
	  axis.ticks = element_blank(),
	  panel.border = element_blank(),
	  strip.background = element_blank(),
	  strip.text = element_text(size = 10, face = "bold"),
	  plot.background = element_rect(color = "white", fill = "white")) +
    labs(fill = NULL)


bcell_prots <- 
    c("IGHM" = "IgM", 
      "FCER2" = "CD23",
      "CD27" = "CD27",
      "CD69" = "CD69", 
      "ITGAX" = "CD11c")

bcell_prot_quant <- bcells@assays$ADT@data |> 
    {function(x) x[bcell_prots, , drop = FALSE]}() |>
    as_tibble(rownames = "gene_id") |>
    pivot_longer(-gene_id, names_to = "barcode", values_to = "prot_exp") |>
    group_by(gene_id) |>
    mutate(prot_exp = clip_expression(prot_exp)) |>
    ungroup() |>
    select(barcode, gene_id, prot_exp) |>
    mutate(gene_id = factor(gene_id, levels = bcell_prots))

bcell_prot_plot_list <- bcell_prot_quant |>
    left_join(umap_df3) |>
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
		  strip.text = element_text(size = 10, face = "bold")) +
	    labs(color = NULL))

bcell_prot_plot <- plot_grid(plotlist = c(bcell_prot_plot_list, list(umap_ki67))) +
    theme(plot.background = element_rect(fill = "white", color = "white"))

ggsave("./plots/bcell_markers_poster.png", bcell_prot_plot, 
       width = 8, height = 4.5, dpi = 600)








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
