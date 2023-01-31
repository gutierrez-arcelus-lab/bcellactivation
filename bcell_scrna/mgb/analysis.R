# single-cell data analysis
library(Seurat)
library(scSHC)

# Data wrangling
library(tidyverse)

# Plotting
library(tidytext)
library(ggridges)
library(RColorBrewer)
library(scico)
library(ggsci)
library(cowplot)

stims <- 
    c(
      "Hashtag6" = "unstim 0h",
      "Hashtag7" = "IL4 24h",
      "Hashtag8" = "IL4 72h",
      "Hashtag9" = "BCR 24h",
      "Hashtag10" = "BCR 72h",
      "Hashtag12" = "TLR7 24h",
      "Hashtag13" = "TLR7 72h",
      "Hashtag14" = "DN2 24h",
      "Hashtag15" = "DN2 72h"
    )

cellranger_dir <- 
    file.path("/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/B_cells_citeseq",
	      "SN0268787/broad/hptmp/curtism/bwh10x/KW10456_Maria",
	      "221101_10X_KW10456_bcl/cellranger-6.1.1/GRCh38/BRI-1984/outs",
	      "filtered_feature_bc_matrix")

features_df <- file.path(cellranger_dir, "features.tsv.gz") |>
    read_tsv(col_names = c("gene_id", "gene_name", "phenotype"))

genes_df <- features_df |>
  filter(phenotype == "Gene Expression") |>
  select(gene_id, gene_name)

mt_genes <- genes_df |>
    filter(grepl("^MT-", gene_name)) |>
    pull(gene_id)

ribo_genes <- genes_df |>
    filter(grepl("^RPS\\d+|^RPL\\d+", gene_name)) |>
    pull(gene_id)

data10x <- Read10X(cellranger_dir, gene.column = 1)

antibody_mtx <- data10x[["Antibody Capture"]] |>
    {function(x) x[!grepl("^Hashtag", rownames(x)), ]}()

rownames(antibody_mtx) <- rownames(antibody_mtx) |>
    sub("_prot$", "", x = _) |>
    gsub("_", ".", x = _)

hashtags_mtx <- data10x[["Antibody Capture"]] |>
    {function(x) x[grepl("^Hashtag", rownames(x)), ]}()

rownames(hashtags_mtx) <- setNames(stims[rownames(hashtags_mtx)], NULL)

stim_colors <- c("grey80", "grey40", "black", brewer.pal(n = 6, "Paired"))
names(stim_colors) <- stims

# Create object
bcells <- CreateSeuratObject(counts = data10x[["Gene Expression"]], project = "MGB")
bcells[["ADT"]] <- CreateAssayObject(counts = antibody_mtx)
bcells[["HTO"]] <- CreateAssayObject(counts = hashtags_mtx)

bcells <- bcells |>
    NormalizeData(normalization.method = "LogNormalize", margin = 1) |>
    NormalizeData(assay = "HTO", normalization.method = "CLR", margin = 1) |>
    NormalizeData(assay = "ADT", normalization.method = "CLR", margin = 2)

bcells[["percent_mt"]] <- PercentageFeatureSet(bcells, features = mt_genes)
bcells[["percent_ribo"]] <- PercentageFeatureSet(bcells, features = ribo_genes)

# Demuxlet
demuxlet_df <- read_tsv("./demuxlet/demuxlet_results.best") |>
    select(barcode = BARCODE, best = BEST) |>
    extract(best, c("status", "sample"), "([^-]+)-(.+)")


# HTO counts
hto <- as_tibble(t(bcells@assays$HTO@counts), rownames = "barcode") |>
    pivot_longer(-barcode, names_to = "stim") |>
    group_by(barcode) |>
    mutate(hto_max = stim[which.max(value)]) |>
    ungroup() |>
    mutate_at(vars(stim, hto_max), ~factor(., levels = stims))

png("./plots/test.png", units = "in", res = 300, width = 6, height = 3)
ggplot(hto, aes(x = log10(value + 1))) +
    geom_density(aes(fill = stim), size = .25, alpha = .5) +
    scale_fill_manual(values = stim_colors) +
    facet_wrap(~hto_max, nrow = 2) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(x = "log10 (counts + 1)")
dev.off()

# QC
bcells <- HTODemux(bcells, assay = "HTO", positive.quantile = 0.99)

bcells_meta <- bcells@meta.data |>
    as_tibble(rownames = "barcode") 

qc_plot <- bcells_meta |>
    select(barcode, n_genes = nFeature_RNA, percent_mt, stim = HTO_maxID) |> 
    ggplot(aes(x = n_genes, y = percent_mt)) +
    geom_point(size = .25, alpha = .25) +
    geom_vline(xintercept = 700, linetype = 2, color = "tomato3") +
    geom_hline(yintercept = 7, linetype = 2, color = "tomato3") +
    facet_wrap(~stim, ncol = 3) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    labs(x = "Number of genes", y = "% reads from Mitochondria")

ggsave("./plots/test.png", qc_plot, width = 6, height = 6)
   

demuxlet_singlets <- demuxlet_df |>
    filter(status == "SNG")

cells_keep <- bcells_meta |>
    filter(nFeature_RNA > 700, percent_mt <= 7) |>
    filter(barcode %in% demuxlet_singlets$barcode)

bcells_pass <- subset(bcells, cells = cells_keep$barcode)

# HTO demux
bcells_pass <- HTODemux(bcells_pass, assay = "HTO", positive.quantile = 0.99)

bcells_pass_demux_meta <- bcells_pass@meta.data |>
    as_tibble(rownames = "barcode") |>
    filter(HTO_classification.global == "Singlet") |>
    inner_join(demuxlet_singlets, by = "barcode")

bcells_singlet <- subset(bcells_pass, cells = bcells_pass_demux_meta$barcode)

## PCA
bcells_singlet <- bcells_singlet |>
    FindVariableFeatures(nfeatures = 2000, selection.method = "vst") |>
    {function(x) ScaleData(x, features = rownames(x))}() |>
    {function(x) RunPCA(x, features = VariableFeatures(x))}()

singlets_meta_data <- bcells_singlet@meta.data |>
    as_tibble(rownames = "barcode") |>
    select(barcode, stim = HTO_maxID, 
	   percent_mt, percent_ribo, n_genes = nFeature_RNA, umi = nCount_RNA) |>
    mutate(stim = factor(stim, levels = stims))

pca_df <- bcells_singlet@reductions$pca@cell.embeddings |>
    as_tibble(rownames = "barcode") |>
    left_join(singlets_meta_data) |>
    select(barcode, stim, n_genes, percent_mt, percent_ribo, PC_1:PC_4)

pca_plot <- ggplot(pca_df, aes(PC_1, PC_2)) +
    geom_point(aes(color = stim), alpha = .5, size = .5) +
    scale_color_manual(values = stim_colors) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 4)))

ggsave("./plots/test.png", pca_plot)

sdev_plot <- tibble(sdev = bcells_singlet@reductions$pca@stdev) |>
    rownames_to_column("pc") |>
    mutate(pc = fct_inorder(pc)) |>
    filter(pc %in% 1:50) |>
    ggplot(aes(pc, sdev)) +
    geom_line(aes(group = 1)) +
    scale_x_discrete(breaks = c(1, seq(5, 25, 5))) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    labs(x = "Principal component", y = "Standard deviation")

loadings_plot <- bcells_singlet@reductions$pca@feature.loadings |>
    as_tibble(rownames = "gene_id") |>
    select(gene_id, PC_1, PC_2) |>
    pivot_longer(-gene_id, names_to = "pc") |>
    mutate(direction = case_when(sign(value) == 1 ~ "+",
				 sign(value) == -1 ~ "-",
				 sign(value) == 0 ~ "0")) |>
    group_by(pc, direction) |>
    top_n(15, abs(value)) |>
    ungroup() |>
    left_join(genes_df) |>
    ggplot(aes(value, reorder_within(gene_name, value, pc))) +
    geom_col(fill = "midnightblue", alpha = .5) +
    scale_y_discrete(labels = function(x) sub("^([^_]+).*$", "\\1", x)) +
    facet_wrap(~pc, scales = "free_y") +
    labs(x = "Loading", y = NULL) +
    theme_bw()

ggsave("./plots/test.png", plot_grid(pca_plot, loadings_plot, nrow = 1), 
       width = 10, height = 5)


bcells_singlet <- bcells_singlet |>
    RunUMAP(dims = 1:30, verbose = FALSE)

umap_df <- as.data.frame(bcells_singlet@reductions$umap@cell.embeddings) |>
    as_tibble(rownames = "barcode") |>
    left_join(singlets_meta_data)

umap_stims <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = stim)) +
    geom_point(size = .5) +
    scale_color_manual(values = stim_colors) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          legend.key.height = unit(.75, "lines")) +
    guides(color = guide_legend(override.aes = list(size = 2)))

mki_gene_df <- genes_df |>
    filter(gene_name == "MKI67") |>
    select(gene_id, gene_name)

mki_gene_quant <- bcells_singlet@assays$RNA@data |> 
    {function(x) x[mki_gene_df$gene_id, ,drop = FALSE]}() |>
    as_tibble(rownames = "gene_id") |>
    left_join(mki_gene_df, by = "gene_id") |>
    pivot_longer(-c("gene_id", "gene_name"), 
                 names_to = "barcode", values_to = "gene_exp")

ki67 <- umap_df |>
    select(barcode, UMAP_1, UMAP_2) |>
    left_join(mki_gene_quant, by = "barcode") |>
    ggplot(aes(UMAP_1, UMAP_2, color = gene_exp, fill = gene_exp)) +
    geom_point(size = .25, shape = 19) +
    scale_color_viridis_c(guide = "none") +
    scale_fill_viridis_c(guide = guide_colorbar(barwidth = .5)) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
	  plot.background = element_rect(fill = "white", color = "white")) +
    labs(title = "MKI67 RNA expression", fill = NULL)

ggsave("./plots/umap.png", plot_grid(umap_stims, ki67, nrow = 1), 
       width = 10, height = 5)


# Marker genes
plot_markers <- function(cluster_df, seurat_obj) {
  
    top_markers <- cluster_df |>
	as_tibble() |>
	left_join(select(features_df, 1:2), by = c("gene" = "gene_id")) |>
	group_by(cluster) |>
	top_n(10, avg_log2FC) |>
	ungroup() |>
	select(cluster, gene_name, gene_id = gene, avg_log2FC)
  
    cluster_cell <- 
	tibble(cluster_cell = Idents(seurat_obj),
	       barcode = names(Idents(seurat_obj))) |>
	mutate(cluster_cell = factor(cluster_cell, levels = levels(top_markers$cluster)))
  
    cell_gene_expr <- seurat_obj@assays$RNA@data |>
        {function(x) x[unique(top_markers$gene_id), ]}() |>
	as.data.frame() |>
	rownames_to_column("gene_id") |>
	as_tibble() |>
	pivot_longer(-gene_id, names_to = "barcode", values_to = "logexpr")
  
    top_marker_expr <- cell_gene_expr |>
	left_join(cluster_cell, by = "barcode") |>
	inner_join(top_markers, by = "gene_id")
  
    top_marker_perc_exp <- top_marker_expr |>
	group_by(cluster = cluster_cell, gene_name) |>
	summarise(prop_expr = mean(logexpr > 0)) |>
	ungroup()
  
    top_marker_avg_exp <- top_marker_expr |>
	filter(logexpr > 0) |>
	group_by(cluster = cluster_cell, gene_name) |>
	summarise(scaled_expr = mean(logexpr)) |>
	ungroup()
  
    top_marker_summary <- left_join(top_marker_perc_exp, top_marker_avg_exp) |>
	left_join(select(top_markers, cluster_top = cluster, gene_name, avg_log2FC)) |>
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

Idents(bcells_singlet) <- "HTO_maxID"

bcells_markers <- 
    FindAllMarkers(bcells_singlet, 
                   only.pos = TRUE,
                   min.pct = 0.1,
                   logfc.threshold = 1) |>
    as_tibble() |>
    mutate(cluster = factor(cluster, levels = stims))

markers_plot <- plot_markers(bcells_markers, bcells_singlet)

ggsave("./plots/markers.png", markers_plot)



# Clustering
scdata <- bcells_singlet@assays$RNA@data

clusters <- scSHC(scdata, 
		  batch = NULL, alpha = 0.05, num_features = 2500,
		  num_PCs = 25, parallel = TRUE, cores = 4)


str(clusters)
