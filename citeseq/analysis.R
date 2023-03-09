# Set large RAM limit for R
unix::rlimit_as(1e12)

# single-cell data analysis
library(Seurat)
library(sclibr)
library(harmony)

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
      "grey50", "grey40",
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

# Import Seurat objects
lib1984 <- read_rds("./data/seurat_1984_qced.rds")
lib1988 <- read_rds("./data/seurat_1988_qced.rds")
lib1990 <- read_rds("./data/seurat_1990_qced.rds")

# Merge datasets and run PCA
bcells <- 
    merge(lib1984, y = c(lib1988, lib1990), add.cell.ids = c("1984", "1988", "1990")) |>
    {function(x) ScaleData(x, features = rownames(x))}() |>
    FindVariableFeatures() |>
    {function(x) RunPCA(x, features = VariableFeatures(x))}()


# Visualize PCA
sdev_plot <- tibble(sdev = bcells@reductions$pca@stdev) |>
    rownames_to_column("pc") |>
    mutate(pc = fct_inorder(pc)) |>
    filter(pc %in% 1:50) |>
    ggplot(aes(pc, sdev)) +
    geom_line(aes(group = 1)) +
    scale_x_discrete(breaks = c(1, seq(5, 50, 5))) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    labs(x = "Principal component", y = "Standard deviation")

ggsave("./plots/pca_sdev.png", sdev_plot, width = 3, height = 2)

# Visualize UMAP
plot_umap <- function(dat, var_name) {

    sc1 <- range(dat$UMAP_1) |>
	{function(x) length(x[1]:x[2])}() |>
	{function(x) x/6}()

    sc2 <- range(dat$UMAP_2) |>
	{function(x) length(x[1]:x[2])}() |>
	{function(x) x/6}()
    
    x1 <- min(dat$UMAP_1)
    x2 <- x1 + sc1
    y1 <- min(dat$UMAP_2)
    y2 <-y1 + sc2 

    ggplot(dat, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(color = {{var_name}}), size = .1, alpha = .25) +
    geom_segment(x = x1, xend = x2, y = y1, yend = y1,
		 arrow = arrow(length = unit(5, "pt")),
		 show.legend = FALSE) +
    geom_segment(x = x1, xend = x1, y = y1, yend = y2,
		 arrow = arrow(length = unit(5, "pt")), 
		 show.legend = FALSE) +
    annotate("text", label = "UMAP 1", size = rel(3),
	     x = x1, y = y1, hjust = 0, vjust = 1.25) + 
    annotate("text", label = "UMAP 2", size = rel(3),
	     x = x1, y = y1, angle = 90, hjust = 0, vjust = -.5) +
    theme_minimal() +
    theme(axis.title = element_blank(),
	  axis.text = element_blank(),
	  panel.grid = element_blank(),
          panel.border = element_blank(),
	  plot.background = element_rect(fill = "white", color = "white")) + 
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))
}

bcells <- bcells |>
    RunUMAP(reduction = "pca", dims = 1:35, reduction.name = "umap.merged")

umap_df <- bcells@reductions$umap.merged@cell.embeddings |>
    as_tibble(rownames = "barcode") |>
    left_join(as_tibble(bcells@meta.data, rownames = "barcode")) |>
    select(barcode, UMAP_1, UMAP_2, dataset = orig.ident, hto) |>
    mutate(hto = factor(hto, levels = stim_order))

umap_stims <- plot_umap(umap_df, hto) +
    scale_color_manual(values = stim_colors) +
    theme(legend.position = c(.8, .25),
          legend.key.height = unit(.75, "lines")) +
    labs(color = "Stim")

umap_batch <- plot_umap(umap_df, dataset) + 
    scale_color_manual(values = c("goldenrod2", "red4", "royalblue")) +
    theme(legend.position = c(.8, .25)) +
    labs(color = "batch")

ggsave("./plots/umap_merge.png", 
       plot_grid(umap_batch, umap_stims, nrow = 1), 
       width = 8, height = 3)


# Run Harmony
bcells <- bcells |>
    RunHarmony("orig.ident") |>
    RunUMAP(reduction = "harmony", dims = 1:35)

# UMAP
umap_df_h <- bcells@reductions$umap@cell.embeddings |>
    as_tibble(rownames = "barcode") |>
    left_join(as_tibble(bcells@meta.data, rownames = "barcode")) |>
    select(barcode, UMAP_1 = umap_1, UMAP_2 = umap_2, dataset = orig.ident, hto) |>
    mutate(hto = factor(hto, levels = stim_order))

umap_stims_h <- plot_umap(umap_df_h, hto) +
    scale_color_manual(values = stim_colors) +
    theme(legend.position = c(.2, .5),
          legend.key.height = unit(.75, "lines")) +
    labs(color = "Stim")

umap_batch_h <- plot_umap(umap_df_h, dataset) + 
    scale_color_manual(values = c("goldenrod2", "red4", "royalblue")) +
    theme(legend.position = c(.2, .5)) +
    labs(color = "batch")

ggsave("./plots/umap_harmony.png", 
       plot_grid(umap_batch_h, umap_stims_h, nrow = 1), 
       width = 8, height = 3)


# Seurat integration with CCA

# Integrate batches
bcell_objects <- list("1984" = lib1984,
		      "1988" = lib1988,
		      "1990" = lib1990)

integration_features <- SelectIntegrationFeatures(object.list = bcell_objects)

anchors <- FindIntegrationAnchors(object.list = bcell_objects,
				  anchor.features = integration_features,
				  assay = c("RNA", "RNA", "RNA"),
				  reduction = "cca",
				  k.anchor = 20,
				  dims = 1:35)

bcells_integrated <- 
    IntegrateData(anchorset = anchors, 
		  dims = 1:35,
		  features.to.integrate = rownames(lib1984@assays$RNA@counts))

write_rds(bcells_integrated, "./data/seurat_mgb_integrated.rds")
#bcells_integrated <- read_rds("./data/seurat_mgb_integrated.rds")

# PCA on integrated data
bcells_integrated <- bcells_integrated |>
    {function(x) ScaleData(x, features = rownames(x))}() |>
    {function(x) RunPCA(x, features = VariableFeatures(x))}() |>
    RunUMAP(dims = 1:35, verbose = FALSE)

# UMAP
umap_df_s <- bcells_integrated@reductions$umap@cell.embeddings |>
    as_tibble(rownames = "barcode") |>
    left_join(as_tibble(bcells_integrated@meta.data, rownames = "barcode")) |>
    select(barcode, UMAP_1, UMAP_2, dataset = orig.ident, hto) |>
    mutate(hto = factor(hto, levels = stim_order))

umap_stims_s <- plot_umap(umap_df_s, hto) +
    scale_color_manual(values = stim_colors) +
    theme(legend.position = c(.75, .25),
          legend.key.height = unit(.75, "lines")) +
    labs(color = "Stim")

umap_batch_s <- plot_umap(umap_df_s, dataset) + 
    scale_color_manual(values = c("goldenrod2", "red4", "royalblue")) +
    theme(legend.position = c(.75, .25)) +
    labs(color = "batch")

ggsave("./plots/umap_seuratIntegration.png", 
       plot_grid(umap_batch_s, umap_stims_s, nrow = 1), 
       width = 8, height = 3)

