
# Set large RAM limit for R
unix::rlimit_as(1e12)

# single-cell data analysis
library(Seurat)
#library(scSHC)


# Data wrangling
library(tidyverse)

# Plotting
library(tidytext)
library(ggridges)
library(RColorBrewer)
library(scico)
library(ggsci)
library(ggthemes)
library(cowplot)

# Function
make_seurat <- function(cellranger_path, project_id) {

    data10x <- Read10X(cellranger_path, gene.column = 1)

    antibody_mtx <- data10x[["Antibody Capture"]] |>
	{function(x) x[!grepl("^Hashtag", rownames(x)), ]}()

    rownames(antibody_mtx) <- sub("_prot$", "", rownames(antibody_mtx))

    hashtags_mtx <- data10x[["Antibody Capture"]] |>
	{function(x) x[grepl("^Hashtag", rownames(x)), ]}()

    rownames(hashtags_mtx) <- setNames(stims[rownames(hashtags_mtx)], NULL)

    bcells <- CreateSeuratObject(counts = data10x[["Gene Expression"]],
				 project = project_id)

    bcells[["ADT"]] <- CreateAssayObject(counts = antibody_mtx)
    bcells[["HTO"]] <- CreateAssayObject(counts = hashtags_mtx)

    bcells <- bcells |>
	NormalizeData(normalization.method = "LogNormalize", margin = 1) |>
	NormalizeData(assay = "HTO", normalization.method = "CLR", margin = 1) |>
	NormalizeData(assay = "ADT", normalization.method = "CLR", margin = 2)

    bcells[["percent_mt"]] <- PercentageFeatureSet(bcells, features = mt_genes)
    bcells[["percent_ribo"]] <- PercentageFeatureSet(bcells, features = ribo_genes)

    bcells
}

# Stims
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

stim_colors <- c("grey80", "grey40", "black", brewer.pal(n = 6, "Paired"))
names(stim_colors) <- stims

# Cell ranger data
libraries <- c("1984", "1988", "1990")

cellranger_dir_1984 <- 
    file.path("/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/B_cells_citeseq",
	      "SN0268787/broad/hptmp/curtism/bwh10x/KW10456_Maria",
	      "221101_10X_KW10456_bcl/cellranger-6.1.1/GRCh38/BRI-1984/outs",
	      "filtered_feature_bc_matrix")

cellranger_dir_1988 <-
    file.path("/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/B_cells_citeseq",
	      "SN0273471/broad/hptmp/sgurajal/bwh10x/KW10598_mgutierrez",
	      "221211_10X_KW10598_bcl/cellranger-6.1.1/GRCh38/BRI-1988_hashing/outs",
	      "filtered_feature_bc_matrix")

cellranger_dir_1990 <-
    file.path("/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/B_cells_citeseq",
	      "SN0273471/broad/hptmp/sgurajal/bwh10x/KW10598_mgutierrez",
	      "221211_10X_KW10598_bcl/cellranger-6.1.1/GRCh38/BRI-1990_hashing/outs",
	      "filtered_feature_bc_matrix")

# Using features from library 1984
# They are the same for all libraries
features_df <- file.path(cellranger_dir_1984, "features.tsv.gz") |>
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

# Seurat objects
bcells_1984 <- make_seurat(cellranger_dir_1984, "1984")
bcells_1988 <- make_seurat(cellranger_dir_1988, "1988")
bcells_1990 <- make_seurat(cellranger_dir_1990, "1990")


# Demuxlet
demuxlet_1984 <- read_tsv("./demuxlet/demuxlet_results.best") |>
    select(barcode = BARCODE, best = BEST) |>
    extract(best, c("status", "sample"), "([^-]+)-(.+)")

demuxlet_1988 <- read_tsv("./demuxlet/demuxlet_1988_results.best") |>
    select(barcode = BARCODE, best = BEST) |>
    extract(best, c("status", "sample"), "([^-]+)-(.+)")

demuxlet_1990 <- read_tsv("./demuxlet/demuxlet_1990_results.best") |>
    select(barcode = BARCODE, best = BEST) |>
    extract(best, c("status", "sample"), "([^-]+)-(.+)")

# HTO counts
plot_hto <- function(seurat_obj) {

    hto <- as_tibble(t(as.matrix(seurat_obj@assays$HTO@counts)), rownames = "barcode") |>
	pivot_longer(-barcode, names_to = "stim") |>
	group_by(barcode) |>
	mutate(hto_max = stim[which.max(value)]) |>
	ungroup() |>
	mutate_at(vars(stim, hto_max), ~factor(., levels = stims))

    ggplot(hto, aes(x = log2(value + 1))) +
	geom_density(aes(fill = stim), linewidth = .25, alpha = .5) +
	scale_fill_manual(values = stim_colors) +
	facet_wrap(~hto_max, nrow = 2) +
	theme_bw() +
	theme(panel.grid = element_blank(),
	      plot.margin = margin(t = 1, unit = "cm")) +
	labs(x = "log2 (counts + 1)")
}


hto_plots <- list(bcells_1984, bcells_1988, bcells_1990) |>  
    map(plot_hto) |>
    plot_grid(plotlist = _, ncol = 1, labels = libraries)

ggsave("./plots/hto.png", hto_plots, width = 6, height = 9)

# QC
plot_qc <- function(seurat_metadata, library_id) {

    seurat_metadata |>
	select(barcode, n_genes = nFeature_RNA, percent_mt, stim = HTO_maxID) |> 
	mutate(stim = factor(stim, levels = stims)) |>
	ggplot(aes(x = n_genes, y = percent_mt)) +
	geom_point(size = .25, alpha = .25) +
	geom_vline(xintercept = 1000, linetype = 2, color = "tomato3") +
	geom_hline(yintercept = 10, linetype = 2, color = "tomato3") +
	facet_wrap(~stim, ncol = 3) +
	theme_bw() +
	theme(panel.grid.minor = element_blank(),
	      plot.margin = margin(t = 1, r = 1, unit = "cm")) +
	labs(x = "Number of genes", 
	     y = "% reads from Mitochondria")
}

bcells_1984 <- HTODemux(bcells_1984, assay = "HTO", positive.quantile = 0.99)
bcells_1988 <- HTODemux(bcells_1988, assay = "HTO", positive.quantile = 0.99)
bcells_1990 <- HTODemux(bcells_1990, assay = "HTO", positive.quantile = 0.99)

bcells_1984_meta <- bcells_1984@meta.data |>
    as_tibble(rownames = "barcode") 

bcells_1988_meta <- bcells_1988@meta.data |>
    as_tibble(rownames = "barcode") 

bcells_1990_meta <- bcells_1990@meta.data |>
    as_tibble(rownames = "barcode") 

qc_plots <- list(bcells_1984_meta, bcells_1988_meta, bcells_1990_meta) |>  
    map(plot_qc) |>
    plot_grid(plotlist = _, ncol = 1, labels = libraries)

ggsave("./plots/qc.png", qc_plots, width = 6, height = 12)
 


# t-SNE
make_tsne_df <- function(x) {
    
    tsne_df <- x@reductions$tsne@cell.embeddings |>
	as_tibble(rownames = "barcode")
    
    meta_df <- as_tibble(x@meta.data, rownames = "barcode") |>
	select(barcode, clas = HTO_classification.global, stim = HTO_maxID)

    left_join(tsne_df, meta_df, by = "barcode")
}


bcells_1984_noneg <- subset(bcells_1984, idents = "Negative", invert = TRUE)
bcells_1988_noneg <- subset(bcells_1988, idents = "Negative", invert = TRUE)
bcells_1990_noneg <- subset(bcells_1990, idents = "Negative", invert = TRUE)

bcells_1984_noneg <- bcells_1984_noneg |>
    {function(x) ScaleData(x, features = rownames(x), verbose = FALSE)}() |>
    FindVariableFeatures() |>
    {function(x) RunPCA(x, features = VariableFeatures(x), verbose = FALSE)}() |>
    {function(x) RunTSNE(x, dims = 1:35, perplexity = 100, nthreads = 4, verbose = TRUE)}()

bcells_1988_noneg <- bcells_1988_noneg |>
    {function(x) ScaleData(x, features = rownames(x), verbose = FALSE)}() |>
    FindVariableFeatures() |>
    {function(x) RunPCA(x, features = VariableFeatures(x), verbose = FALSE)}() |>
    {function(x) RunTSNE(x, dims = 1:35, perplexity = 100, nthreads = 4, verbose = TRUE)}()

bcells_1990_noneg <- bcells_1990_noneg |>
    {function(x) ScaleData(x, features = rownames(x), verbose = FALSE)}() |>
    FindVariableFeatures() |>
    {function(x) RunPCA(x, features = VariableFeatures(x), verbose = FALSE)}() |>
    {function(x) RunTSNE(x, dims = 1:35, perplexity = 100, nthreads = 4, verbose = TRUE)}()

tsne_dat <- bind_rows(
    "1984" = make_tsne_df(bcells_1984_noneg),
    "1988" = make_tsne_df(bcells_1988_noneg),
    "1990" = make_tsne_df(bcells_1990_noneg),
    .id = "set"
)


demuxlet_df <- 
    bind_rows("1984" = demuxlet_1984,
	      "1988" = demuxlet_1988,
	      "1990" = demuxlet_1990, 
	      .id = "set") |>
    select(set, barcode, status)

stim_colors_v2 <- c("goldenrod1", "goldenrod3", "goldenrod4", brewer.pal(n = 6, "Paired"))
names(stim_colors_v2) <- stims
stim_colors_v2 <- c("doublet" = "black", stim_colors_v2)

doublet_df <- tsne_dat |>
    left_join(demuxlet_df)

doublet_hto <- doublet_df |>
    mutate(doublet_lab = case_when(clas == "Doublet" & status == "SNG" ~ "doublet",
				   TRUE ~ stim))

doublet_demux <- doublet_df |>
    mutate(doublet_lab = case_when(clas != "Doublet" & status != "SNG" ~ "doublet",
				   TRUE ~ stim))

doublet_both <- doublet_df |>
    mutate(doublet_lab = case_when(clas == "Doublet" & status != "SNG" ~ "doublet",
				   TRUE ~ stim))


doublet_plot_df <- 
    bind_rows("Doublets for HTO only" = doublet_hto,
	      "Doublets for demuxlet only" = doublet_demux,
	      "Doublets for both" = doublet_both,
	      .id = "rowid") |>
    mutate(rowid = factor(rowid, levels = c("Doublets for HTO only",
					    "Doublets for demuxlet only",
					    "Doublets for both")),
	   doublet_lab = factor(doublet_lab, levels = names(stim_colors_v2)))
    
tsne_tmp <- ggplot() +
    geom_point(data = doublet_plot_df, 
	       aes(tSNE_1, tSNE_2, color = doublet_lab), 
	       size = 4, alpha = 1) +
    scale_color_manual(values = stim_colors_v2) +
    labs(color = "Condition")

tsne_legend <- get_legend(tsne_tmp)

tsne_plot <- ggplot() +
    geom_point(data = filter(doublet_plot_df, !grepl("doublet", doublet_lab)),
	       aes(tSNE_1, tSNE_2, color = doublet_lab), size = .2, alpha = .1) +
    geom_point(data = filter(doublet_plot_df, grepl("doublet", doublet_lab)),
	       aes(tSNE_1, tSNE_2, color = doublet_lab), size = .2, alpha = .5) +
    scale_color_manual(values = stim_colors_v2) +
    facet_grid(rowid~set) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
	  axis.title = element_text(size = 8),
	  plot.background = element_rect(fill = "white", color = "white"),
	  legend.position = "none")


ggsave("./plots/tsne.png", 
       plot_grid(tsne_plot, tsne_legend, rel_widths = c(1, .2)) + 
	   theme(plot.background = element_rect(fill = "white", color = "white")), 
       width = 8, height = 6)

doublet_counts_plot <- 
    doublet_df |>
    mutate(lab = case_when(clas == "Doublet" & status != "SNG" ~ "Doublet for both",
			   clas == "Doublet" & status == "SNG" ~ "Doublet for HTO",
			   clas == "Singlet" & status != "SNG" ~ "Doublet for demuxlet",
			   TRUE ~ "Singlet")) |>
    count(set, lab) |>
    group_by(set) |>
    mutate(prop = n/sum(n)) |>
    ungroup() |>
    mutate(lab = factor(lab, levels = c("Singlet", "Doublet for HTO", "Doublet for demuxlet", "Doublet for both"))) |>
    ggplot(aes(lab, prop)) +
	geom_col() +
	facet_wrap(~set, nrow = 1) +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
	labs(x = NULL, y = "Proportion of droplets")

ggsave("./plots/doublet_counts.png", doublet_counts_plot, width = 7, height = 3.5)

# Plot admix
library(sclibr)

admix_1984_s <- subset(bcells_1984_noneg, idents = "Doublet", invert = TRUE) |>
    plot_admix(stim_colors) +
    labs(title = "1984: Singlets")

admix_1984_d <- subset(bcells_1984_noneg, idents = "Doublet") |>
    plot_admix(stim_colors) +
    labs(title = "1984: Doublets")

ggsave("./plots/admix_1984.png", 
       plot_grid(admix_1984_s, admix_1984_d, ncol = 1),
       width = 8, height = 4)

admix_1988_s <- subset(bcells_1988_noneg, idents = "Doublet", invert = TRUE) |>
    plot_admix(stim_colors) +
    labs(title = "1988: Singlets")

admix_1988_d <- subset(bcells_1988_noneg, idents = "Doublet") |>
    plot_admix(stim_colors) +
    labs(title = "1988: Doublets")

ggsave("./plots/admix_1988.png", 
       plot_grid(admix_1988_s, admix_1988_d, ncol = 1),
       width = 8, height = 4)

admix_1990_s <- subset(bcells_1990_noneg, idents = "Doublet", invert = TRUE) |>
    plot_admix(stim_colors) +
    labs(title = "1990: Singlets")

admix_1990_d <- subset(bcells_1990_noneg, idents = "Doublet") |>
    plot_admix(stim_colors) +
    labs(title = "1990: Doublets")

ggsave("./plots/admix_1990.png", 
       plot_grid(admix_1990_s, admix_1990_d, ncol = 1),
       width = 8, height = 4)







# Filter good cells
filter_cells <- function(seurat_obj, n_genes = 1000, p_mt = 10, demuxlet_df) {

    demuxlet_singlets <- demuxlet_df |>
	filter(status == "SNG") |>
	group_by(sample) |>
	filter(n() > 10) |>
	ungroup() |>
	pull(barcode)

    cells_keep <- seurat_obj@meta.data |>
	as_tibble(rownames = "barcode") |>
	filter(nFeature_RNA >= n_genes, 
	       percent_mt <= p_mt,
	       HTO_classification.global == "Singlet",
	       barcode %in% demuxlet_singlets)

    seurat_obj <- subset(seurat_obj, cells = cells_keep$barcode)

    seurat_obj
}

bcells_1984_sng <- filter_cells(bcells_1984, demuxlet_df = demuxlet_1984)
bcells_1988_sng <- filter_cells(bcells_1988, demuxlet_df = demuxlet_1988)
bcells_1990_sng <- filter_cells(bcells_1990, demuxlet_df = demuxlet_1990)


# Plot cells per donor
cells_before_qc <- 
    bind_rows("1984" = demuxlet_1984, 
	      "1988" = demuxlet_1988, 
	      "1990" = demuxlet_1990,
	      .id = "library_id") |>
    filter(status == "SNG") |>
    select(-status)

cells_after_qc <- 
    map_df(setNames(libraries, libraries),
	    function(x) {
		seurat_obj <- get(sprintf("bcells_%s_sng", x)) 
		demux_df <- get(sprintf("demuxlet_%s", x))
		
		seurat_obj@meta.data |>
		as_tibble(rownames = "barcode") |>
		select(barcode) |>
		left_join(demux_df) |>
		select(-status)
	    }, 
	    .id = "library_id")


n_cells_plot <- 
    bind_rows("Before QC" = cells_before_qc, "After QC" = cells_after_qc, .id = "qc") |>
    count(library_id, qc, sample) |>
    extract(sample, "sample", ".+-(\\d+)") |>
    mutate(qc = factor(qc, levels = c("Before QC", "After QC"))) |>
    ggplot(aes(x = sample, y = n)) +
	geom_col(fill = "midnightblue", alpha = .5) +
	scale_y_continuous(breaks = seq(from = 0, to = 16000, by = 2000)) +
	facet_grid(qc ~ library_id, scales = "free_x") +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
	      panel.grid.minor.y = element_blank()) +
	labs(x = NULL, y = "Total cells")

ggsave("./plots/ncells.png", n_cells_plot, width = 4, height = 3)


# Clustering
#scdata <- bcells_singlet@assays$RNA@counts
#
#clusters <- scSHC(scdata, 
#		  batch = NULL, alpha = 0.05, num_features = 2000,
#		  num_PCs = 30, parallel = TRUE, cores = 4)
#
#cluster_df <- enframe(clusters[[1]], "barcode", "cluster")
#
#metadata_update <- bcells_singlet@meta.data |> 
#    as_tibble(rownames = "barcode") |>
#    left_join(cluster_df, by = "barcode") |>
#    column_to_rownames("barcode") |>
#    as.data.frame()
#
#bcells_singlet <- AddMetaData(bcells_singlet, metadata = metadata_update)


# Integrate batches
bcell_objects <- list("1984" = bcells_1984_sng,
		      "1988" = bcells_1988_sng,
		      "1990" = bcells_1990_sng)

integration_features <- SelectIntegrationFeatures(object.list = bcell_objects)

anchors <- FindIntegrationAnchors(object.list = bcell_objects,
				  anchor.features = integration_features,
				  assay = c("RNA", "RNA", "RNA"),
				  reduction = "cca",
				  k.anchor = 20,
				  dims = 1:30)

bcells_integrated <- 
    IntegrateData(anchorset = anchors, 
		  dims = 1:30,
		  features.to.integrate = rownames(bcells_1984_sng@assays$RNA@counts))

# PCA on integrated data
bcells_integrated <- bcells_integrated |>
    {function(x) ScaleData(x, features = rownames(x))}() |>
    {function(x) RunPCA(x, features = VariableFeatures(x))}()

write_rds(bcells_integrated, "./bcells_integrated_seurat.rds")
#bcells_integrated <- read_rds("./bcells_integrated_seurat.rds")

# Plot PCA
integrated_meta <- bcells_integrated@meta.data |>
    as_tibble(rownames = "barcode") |>
    select(barcode, library_id = orig.ident, stim = HTO_maxID, 
	   percent_mt, percent_ribo, n_genes = nFeature_RNA, umi = nCount_RNA) |>
    mutate(stim = factor(stim, levels = stims))

pca_df <- bcells_integrated@reductions$pca@cell.embeddings |>
    as_tibble(rownames = "barcode") |>
    select(barcode, PC_1:PC_6) |>
    left_join(integrated_meta) |>
    extract(barcode, "barcode", "(.+)_[123]") |>
    left_join(cells_before_qc) |>
    extract(sample, "sample_id", "\\d+_[A-Z0-9]+-(\\d+)") |>
    select(barcode, library_id, sample_id, 
	   stim, n_genes, percent_mt, percent_ribo, PC_1:PC_6)


pca_plot1 <- ggplot(pca_df, aes(PC_1, PC_2)) +
    geom_point(aes(color = library_id), alpha = .5, size = .5) +
    scale_color_manual(values = c("goldenrod3", "red", "midnightblue")) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 4)))

pca_plot2 <- ggplot(pca_df, aes(PC_1, PC_2)) +
    geom_point(aes(color = sample_id), alpha = .5, size = .5) +
    scale_color_aaas() +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 4)))

pca_plot3 <- ggplot(pca_df, aes(PC_1, PC_2)) +
    geom_point(aes(color = stim), alpha = .5, size = .5) +
    scale_color_manual(values = stim_colors) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 4)))

sdev_plot <- tibble(sdev = bcells_integrated@reductions$pca@stdev) |>
    rownames_to_column("pc") |>
    mutate(pc = fct_inorder(pc)) |>
    filter(pc %in% 1:50) |>
    ggplot(aes(pc, sdev)) +
    geom_line(aes(group = 1)) +
    scale_x_discrete(breaks = c(1, seq(5, 50, 5))) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    labs(x = "Principal component", y = "Standard deviation")

ggsave("./plots/pca_sdev.png", sdev_plot)

loadings_plot <- bcells_integrated@reductions$pca@feature.loadings |>
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
    theme_bw() +
    theme(axis.text.y = element_text(size = 8))

ggsave("./plots/pca.png", 
       plot_grid(pca_plot1, pca_plot2, pca_plot3, loadings_plot, nrow = 2),
       width = 10)



# UMAP
bcells_integrated <- bcells_integrated |>
    RunUMAP(dims = 1:35, verbose = FALSE)

umap_df <- as.data.frame(bcells_integrated@reductions$umap@cell.embeddings) |>
    as_tibble(rownames = "barcode") |>
    left_join(integrated_meta) |>
    extract(barcode, "barcode", "(.+)_[123]") |>
    left_join(cells_before_qc) |>
    extract(sample, "sample_id", "\\d+_[A-Z0-9]+-(\\d+)")

umap_stims <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = stim)) +
    geom_point(size = .25, alpha = .5) +
    scale_color_manual(values = stim_colors) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          legend.key.height = unit(.75, "lines")) +
    guides(color = guide_legend(override.aes = list(size = 2)))

umap_batch <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = library_id)) +
    geom_point(size = .25, alpha = .5) +
    scale_color_manual(values = c("goldenrod3", "red", "midnightblue")) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          legend.key.height = unit(.75, "lines")) +
    guides(color = guide_legend(override.aes = list(size = 2))) +
    labs(color = "batch")

umap_donor <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = sample_id)) +
    geom_point(size = .25, alpha = .5) +
    scale_color_manual(values = c("black", "goldenrod3", "cornflowerblue", "green", "magenta")) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          legend.key.height = unit(.75, "lines")) +
    guides(color = guide_legend(override.aes = list(size = 2))) +
    labs(color = "donor")

mki_gene_df <- genes_df |>
    filter(gene_name == "MKI67") |>
    select(gene_id, gene_name)

mki_gene_quant <- bcells_integrated@assays$integrated@data |> 
    {function(x) x[mki_gene_df$gene_id, ,drop = FALSE]}() |>
    as_tibble(rownames = "gene_id") |>
    left_join(mki_gene_df, by = "gene_id") |>
    pivot_longer(-c("gene_id", "gene_name"), 
                 names_to = "barcode", values_to = "gene_exp")

ki67 <- as.data.frame(bcells_integrated@reductions$umap@cell.embeddings) |>
    as_tibble(rownames = "barcode") |> 
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
    labs(fill = "KI67\nexpression")

ggsave("./plots/umap.png", 
       plot_grid(umap_stims, umap_batch, umap_donor, ki67, nrow = 2), 
       width = 10, height = 8)




# Marker genes recommended for integrated data
DefaultAssay(bcells_integrated) <- "RNA"
Idents(bcells_integrated) <- "HTO_maxID"


bcr24_markers <- 
    FindConservedMarkers(bcells_integrated, 
			 ident.1 = "BCR 24h",
			 ident.2 = "IL4 24h",
			 grouping.var = "orig.ident", 
			 verbose = FALSE)

dn72_markers <- 
    FindConservedMarkers(bcells_integrated, 
			 ident.1 = "DN2 72h",
			 grouping.var = "orig.ident", 
			 verbose = FALSE)

bcr24_markers |>
    as_tibble(rownames = "gene_id") |>
    left_join(genes_df) |>
    print(width = Inf)
    
dn72_markers |>
    as_tibble(rownames = "gene_id") |>
    left_join(genes_df) |>
    print(width = Inf)


sle_genes <- read_tsv("../reported_genes.tsv")

sle_genes_df <- filter(genes_df, gene_name %in% sle_genes$gene)

b_cell_df <- genes_df |>
    filter(gene_name %in% c("CD27", "TNFRSF13B", "TNFRSF17", "TCL1A", "FCER2", 
			    "IL4R", "JCHAIN", "CD69", "CD83", "ITGAX", "TBX21", "MZB1"))

bcell_gene_quant <- bcells_integrated@assays$integrated@data |> 
    {function(x) x[rownames(x) %in% b_cell_df$gene_id, , drop = FALSE]}() |>
    as_tibble(rownames = "gene_id") |>
    left_join(b_cell_df, by = "gene_id") |>
    pivot_longer(-c("gene_id", "gene_name"), 
                 names_to = "barcode", values_to = "gene_exp")

sle_gene_quant <- bcells_integrated@assays$integrated@data |> 
    {function(x) x[rownames(x) %in% sle_genes_df$gene_id, , drop = FALSE]}() |>
    as_tibble(rownames = "gene_id") |>
    left_join(sle_genes_df, by = "gene_id") |>
    pivot_longer(-c("gene_id", "gene_name"), 
                 names_to = "barcode", values_to = "gene_exp")

sle_gene_umap_df <- as.data.frame(bcells_integrated@reductions$umap@cell.embeddings) |>
    as_tibble(rownames = "barcode") |> 
    right_join(sle_gene_quant, by = "barcode") 

bcell_umap_df <- as.data.frame(bcells_integrated@reductions$umap@cell.embeddings) |>
    as_tibble(rownames = "barcode") |> 
    right_join(bcell_gene_quant, by = "barcode") 


sle_expressed <- sle_gene_umap_df |> 
    group_by(gene_id) |> 
    mutate(avg = mean(gene_exp)) |>
    filter(avg > .1) |>
    ungroup() |>
    distinct(gene_id, gene_name, avg) |>
    arrange(desc(avg))


test_umap <- sle_gene_umap_df |>
    filter(gene_name %in% c("SOCS1", "PTPRC", "IRF8", "BANK1", "NR4A3", "BLK")) |>
    group_split(gene_name) |>
    map(function(x) ggplot(data = x, aes(UMAP_1, UMAP_2, color = gene_exp, fill = gene_exp)) +
    geom_point(size = .25, shape = 19) +
    scale_color_viridis_c(option = "inferno", guide = "none") +
    scale_fill_viridis_c(option = "inferno", guide = guide_colorbar(barwidth = .5)) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
	  plot.background = element_rect(fill = "white", color = "white")) +
    labs(title = unique(x$gene_name), fill = "RNA")) |>
    plot_grid(plotlist = _)

ggsave("./plots/test_umap.png", test_umap, width = 10, height = 7) 


test_umap <- bcell_umap_df |>
    group_split(gene_name) |>
    map(function(x) ggplot(data = x, aes(UMAP_1, UMAP_2, color = gene_exp, fill = gene_exp)) +
    geom_point(size = .25, shape = 19) +
    scale_color_viridis_c(option = "inferno", guide = "none") +
    scale_fill_viridis_c(option = "inferno", guide = guide_colorbar(barwidth = .5)) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
	  plot.background = element_rect(fill = "white", color = "white")) +
    labs(title = unique(x$gene_name), fill = "RNA")) |>
    plot_grid(plotlist = _) +
    theme(plot.background = element_rect(fill = "white", color = "white"))

ggsave("./plots/test_umap.png", test_umap, width = 12, height = 9) 



# default seurat clustering
DefaultAssay(bcells_integrated) <- "integrated"

bcells_integrated <- bcells_integrated |>
    FindNeighbors(dims = 1:35) |>
    FindClusters(resolution = 0.15)

cluster_df <- bcells_integrated@meta.data |>
    as_tibble(rownames = "barcode") |>
    select(barcode, cluster = integrated_snn_res.0.15)

umap_only <- as.data.frame(bcells_integrated@reductions$umap@cell.embeddings) |>
    as_tibble(rownames = "barcode")  

cluster_labs <- umap_only |> 
    left_join(cluster_df, by = "barcode") |>
    group_by(cluster) |>
    summarise_at(vars(UMAP_1, UMAP_2), median) |>
    ungroup()

umap_cluster <- umap_only |> 
    left_join(cluster_df, by = "barcode") |>
    ggplot(aes(UMAP_1, UMAP_2)) +
    geom_point(aes(color = cluster), size = .25) +
    geom_label(data = cluster_labs, aes(label = cluster), alpha = .5) +
    scale_color_manual(values = c("black", "grey", pal_npg()(10))) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
	  plot.background = element_rect(fill = "white", color = "white"),
	  legend.position = "none")

ggsave("./plots/umap_cluster.png", umap_cluster, width = 5, height = 5) 






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

Idents(bcells_integrated) <- "integrated_snn_res.0.15"

cluster_markers <- 
    FindAllMarkers(bcells_integrated, 
                   only.pos = TRUE,
                   min.pct = 0.1,
                   logfc.threshold = .5) |>
    as_tibble()

cluster_markers_plot <- plot_markers(cluster_markers, bcells_integrated)

ggsave("./plots/cluster_markers.png", cluster_markers_plot)


# ADT integration
# Integrate batches
bcell_objects <- list("1984" = bcells_1984_sng,
		      "1988" = bcells_1988_sng,
		      "1990" = bcells_1990_sng)

DefaultAssay(bcells_1984_sng) <- "ADT"
DefaultAssay(bcells_1988_sng) <- "ADT"
DefaultAssay(bcells_1990_sng) <- "ADT"

seurat_adt_list <- 
  list("1984" = bcells_1984_sng, "1988" = bcells_1988_sng, "1900" = bcells_1990_sng) |>
  map(function(x) {
    x <- FindVariableFeatures(x, loess.span = .5)
    x <- ScaleData(x, features = rownames(x), verbose = FALSE)
    x <- RunPCA(x, features = VariableFeatures(x), assay = "ADT",
                approx = FALSE, verbose = FALSE)
  })

adt_anchors <- 
  FindIntegrationAnchors(object.list = seurat_adt_list, 
			 assay = c("ADT", "ADT", "ADT"),
			 reduction = "rpca",
                         k.anchor = 20)

adt_integrated <- 
  IntegrateData(anchorset = adt_anchors, 
                features.to.integrate = rownames(bcells_1984_sng))

adt_df <- adt_integrated@assays$integrated@data |>
    as_tibble(rownames = "ab") |>
    pivot_longer(-ab, names_to = "barcode", values_to = "ab_level") |>
    mutate(ab = factor(ab, levels = str_sort(unique(ab), numeric = TRUE))) |>
    left_join(umap_only)

# shrink outliers
adt_df_sh <- adt_df |>
    group_by(ab) |>
    mutate(q01 = quantile(ab_level, 0.01),
	   q99 = quantile(ab_level, 0.99),
	   ab_level = case_when(ab_level < q01 ~ q01,
				ab_level > q99 ~ q99,
				TRUE ~ ab_level)) |>
    ungroup()

bcell_prots <- 
    c("IGHD" = "IgD", 
      "IGHM" = "IgM", 
      "CR2" = "CD21", 
      "FCER2" = "CD23",
      "CD27" = "CD27",
      "CD69" = "CD69", 
      "CD86" = "CD86",
      "ITGAX" = "CD11c", 
      "CXCR5" = "CD185-or-CXCR5")


bcell_adt_plot <- adt_df_sh |>
    filter(ab %in% bcell_prots) |>
    mutate(ab = as.character(ab)) |>
    {function(x) split(x, x$ab)}() |>
    map(~ggplot(data = ., aes(UMAP_1, UMAP_2, color = ab_level)) +
	    geom_point(size = .2) +
	    scale_color_scico(palette = "lajolla",
			      labels = function(x) str_pad(x, 3),
			      guide = guide_colorbar(barwidth = .5,
						     barheight = 3)) +
	    facet_wrap(~ab) +
	    theme_bw() +
	    theme(panel.grid = element_blank(),
		  panel.border = element_blank(),
		  axis.title = element_blank(),
		  axis.text = element_blank(),
		  strip.background = element_blank(),
		  axis.line = element_blank(),
		  axis.ticks = element_blank()) +
	    labs(color = NULL)) |>
    plot_grid(plotlist = _)

ggsave("./plots/adt.png", bcell_adt_plot, width = 7, height = 6)










# ADT
adt_df <- bcells_singlet@assays$ADT@data |>
    as_tibble(rownames = "ab") |>
    pivot_longer(-ab, names_to = "barcode", values_to = "ab_level") |>
    mutate(ab = factor(ab, levels = str_sort(unique(ab), numeric = TRUE)))

# shrink outliers
adt_df_sh <- adt_df |>
    group_by(ab) |>
    mutate(q01 = quantile(ab_level, 0.01),
	   q99 = quantile(ab_level, 0.99),
	   ab_level = case_when(ab_level < q01 ~ q01,
				ab_level > q99 ~ q99,
				TRUE ~ ab_level)) |>
    ungroup()

bcell_prots_plot_list <- umap_df |>
    select(barcode, UMAP_1, UMAP_2) |>
    left_join(adt_df_sh, by = "barcode", multiple = "all") |>
    {function(x) split(x, x$ab)}() |>
    map(~ggplot(data = ., aes(UMAP_1, UMAP_2, color = ab_level)) +
	    geom_point(size = .2) +
	    scale_color_scico(palette = "lajolla",
			      labels = function(x) str_pad(x, 3),
			      guide = guide_colorbar(barwidth = .5,
						     barheight = 3)) +
	    facet_wrap(~ab) +
	    theme_bw() +
	    theme(panel.grid = element_blank(),
		  panel.border = element_blank(),
		  axis.title = element_blank(),
		  axis.text = element_blank(),
		  strip.background = element_blank(),
		  axis.line = element_blank(),
		  axis.ticks = element_blank()) +
	    labs(color = NULL))

adt_plot <- plot_grid(plotlist = bcell_prots_plot_list, ncol = 4)

ggsave("./plots/adt.png", width = 10, height = 40)




# cell types 

library(MCPcounter)

mcp_scores <- 
    MCPcounter.estimate(expression = bcells_singlet@assays$RNA@data,
		      featuresType = "ENSEMBL_ID") |>
    t() |>
    as.data.frame() |>
    rownames_to_column("barcode") |>
    as_tibble() |>
    pivot_longer(-barcode, names_to = "cell_type")

umap_celltype_list <- umap_df |>
    left_join(mcp_scores, by = "barcode") |>
    {function(x) split(x, x$cell_type)}() |>
    map(~ggplot(.) +
    geom_point(aes(UMAP_1, UMAP_2, color = value, fill = value), size = .1, shape = 19) +
    scale_color_viridis_c(guide = "none") +
    scale_fill_viridis_c("", guide = guide_colorbar(barwidth = .5)) +
    facet_wrap(~cell_type, scales = "free", ncol = 3) +
    theme_bw() +
    theme(panel.grid = element_blank(),
	  panel.border = element_blank(),
	  axis.line = element_blank(),
	  axis.ticks = element_blank(),
	  legend.key.height = unit(.75, "lines"),
	  axis.title = element_blank(),
	  axis.text = element_blank()))

mcp_plot <- plot_grid(plotlist = umap_celltype_list, ncol = 3)

ggsave("./plots/mcp.png", mcp_plot)
