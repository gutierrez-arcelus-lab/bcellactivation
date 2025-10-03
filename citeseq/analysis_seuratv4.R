# RAM
unix::rlimit_as(1e12)

# Packages
## single-cell data analysis
library(Seurat)
library(demuxmix)
library(harmony)

## Data wrangling
library(tidyverse)

# Colors
stim_colors <- 
    "../figure_colors.txt" |>
    read_tsv(col_names = c("stim", "time", "color")) |>
    mutate(stim = glue::glue("{stim} {time}h")) |>
    select(stim, color) |>
    deframe()

# Processing function
make_seurat <- function(cellranger_path, project_id, hto_names = NULL, mito_ids, ribo_ids) {

    data10x <- Read10X(cellranger_path, gene.column = 1)

    antibody_mtx <- data10x[["Antibody Capture"]] |>
        {function(x) x[!grepl("^Hashtag", rownames(x)), ]}()

    if (!is.null(hto_names)) { 

	hashtags_mtx <- data10x[["Antibody Capture"]] |>
	    {function(x) x[grepl("^Hashtag", rownames(x)), ]}()

	rownames(hashtags_mtx) <- setNames(hto_names[rownames(hashtags_mtx)], NULL)
    }

    seuratobj <- 
	CreateSeuratObject(counts = data10x[["Gene Expression"]], project = project_id) |>
        NormalizeData(assay = "RNA", normalization.method = "LogNormalize")

    if ( nrow(antibody_mtx) > 0 ) {
	
	rownames(antibody_mtx) <- sub("_prot$", "", rownames(antibody_mtx))
	
	seuratobj[["ADT"]] <- CreateAssayObject(counts = antibody_mtx)
    
	seuratobj <- NormalizeData(seuratobj, assay = "ADT", normalization.method = "CLR", margin = 2)
    }
    
    if (!is.null(hto_names)) {
	
	seuratobj[["HTO"]] <- CreateAssayObject(counts = hashtags_mtx)

	seuratobj <- NormalizeData(seuratobj, assay = "HTO", normalization.method = "CLR", margin = 2)
    }

    seuratobj[["percent_mt"]] <- PercentageFeatureSet(seuratobj, features = mito_ids)
    seuratobj[["percent_ribo"]] <- PercentageFeatureSet(seuratobj, features = ribo_ids)

    seuratobj
}

run_demuxmix <- function(x) {
   
    hto <- as.matrix(x@assays$HTO@counts)

    rna <- x@meta.data |> 
	select(nFeature_RNA) |>
	as_tibble(rownames = "barcode") |>
	deframe()

    rna <- rna[colnames(hto)]
    
    dmm <- demuxmix(hto, rna = rna, maxIter = 500, tol = 1e-4)
    
    classes <- 
	dmmClassify(dmm) |>
	as_tibble(rownames = "barcode") |>
	select(barcode, dmm_hto_call = HTO, dmm_prob = Prob, dmm_type = Type)
    
    x@meta.data <- 
	x@meta.data |>
	as_tibble(rownames = "barcode") |>
	left_join(classes, join_by(barcode)) |>
	column_to_rownames("barcode")

    x
}

# Directories
labshr <- "/lab-share/IM-Gutierrez-e2/Public/Lab_datasets"

lib1984_dir <- 
    file.path(labshr, "B_cells_citeseq",
	      "SN0268787/broad/hptmp/curtism/bwh10x/KW10456_Maria",
	      "221101_10X_KW10456_bcl/cellranger-6.1.1/GRCh38/BRI-1984/outs",
	      "filtered_feature_bc_matrix")

lib1988_dir <-
    file.path(labshr, "B_cells_citeseq",
	      "SN0273471/broad/hptmp/sgurajal/bwh10x/KW10598_mgutierrez",
	      "221211_10X_KW10598_bcl/cellranger-6.1.1/GRCh38/BRI-1988_hashing/outs",
	      "filtered_feature_bc_matrix")

lib1990_dir <-
    file.path(labshr, "B_cells_citeseq",
	      "SN0273471/broad/hptmp/sgurajal/bwh10x/KW10598_mgutierrez",
	      "221211_10X_KW10598_bcl/cellranger-6.1.1/GRCh38/BRI-1990_hashing/outs",
	      "filtered_feature_bc_matrix")

# Feature IDs
features <- file.path(lib1984_dir, "features.tsv.gz") |>
    read_tsv(col_names = c("gene_id", "gene_name", "phenotype"))

# Mitochondrial and Ribosomal protein genes
# Gene IDs are the same across CITE-seq runs, so we'll use pilot 1.
mt_genes <- features |>
    filter(grepl("^MT-", gene_name)) |>
    pull(gene_id)

ribo_genes <- features |>
    filter(grepl("^RPS\\d+|^RPL\\d+", gene_name)) |>
    pull(gene_id)

# Hashtag IDs
#pilot1_stims <- 
#    c("Hashtag1" = "BCR 72h",
#      "Hashtag2" = "TLR7 72h",
#      "Hashtag3" = "BCR 24h", 
#      "Hashtag4" = "TLR7 24h",
#      "Hashtag5" = "Unstim 24h",
#      "Hashtag6" = "Unstim 0h")
#
#pilot2_stims <- 
#    c("Hashtag1" = "Unstim 0h", 
#      "Hashtag2" = "IL4 24h",
#      "Hashtag3" = "BCR 24h",
#      "Hashtag4" = "BCR-TLR7 24h",
#      "Hashtag5" = "TLR7 24h", 
#      "Hashtag6" = "CD40L 24h",
#      "Hashtag7" = "TLR9 24h",
#      "Hashtag8" = "DN2 24h",
#      "Hashtag9" = "BCR 72h",
#      "Hashtag10" = "BCR-TLR7 72h",
#      "Hashtag12" = "TLR7 72h",
#      "Hashtag13" = "CD40L 72h",
#      "Hashtag14" = "DN2 72h",
#      "Hashtag15" = "TLR9 72h")

mgb_stims <- 
    c("Hashtag6" = "Unstim 0h",
      "Hashtag7" = "IL4 24h",
      "Hashtag8" = "IL4 72h",
      "Hashtag9" = "BCR 24h",
      "Hashtag10" = "BCR 72h",
      "Hashtag12" = "TLR7 24h",
      "Hashtag13" = "TLR7 72h",
      "Hashtag14" = "DN2 24h",
      "Hashtag15" = "DN2 72h")

# Seurat objects
# Demultiplex by HTO using demuxmix
lib1984_obj <- 
    make_seurat(lib1984_dir, 
		project_id = "1984", 
		hto_names = mgb_stims,
		mito_ids = mt_genes, ribo_ids = ribo_genes) |>
    run_demuxmix()

lib1988_obj <- 
    make_seurat(lib1988_dir, 
		project_id = "1988", 
		hto_names = mgb_stims,
		mito_ids = mt_genes, ribo_ids = ribo_genes) |>
    run_demuxmix()

lib1990_obj <- 
    make_seurat(lib1990_dir, 
		project_id = "1990", 
		hto_names = mgb_stims,
		mito_ids = mt_genes, ribo_ids = ribo_genes) |>
    run_demuxmix()


# Filter out droplets with too few genes and %mt > 10
meta_df <- 
    list(lib1984_obj, lib1988_obj, lib1990_obj) |>
    map_dfr(function(x) x@meta.data |> 
	    as_tibble(rownames = "barcode") |>
	    select(barcode, orig.ident, n_genes = nFeature_RNA, n_umi = nCount_RNA, percent_mt))

meta_good <- meta_df |>
    filter(n_genes >= 500, percent_mt <= 10)

good_cells <- split(meta_good, meta_good$orig.ident) |>
    map("barcode")

lib1984_obj <- subset(lib1984_obj, cells = good_cells[["1984"]])
lib1988_obj <- subset(lib1988_obj, cells = good_cells[["1988"]])
lib1990_obj <- subset(lib1990_obj, cells = good_cells[["1990"]])

## Export data for Scrublet
#DropletUtils::write10xCounts(x = lib1984_obj@assays$RNA@counts, 
#			     path = "./data/lib1984_matrix", 
#			     overwrite = TRUE)
#
#DropletUtils::write10xCounts(x = lib1988_obj@assays$RNA@counts, 
#			     path = "./data/lib1988_matrix",
#			     overwrite = TRUE)
#
#DropletUtils::write10xCounts(x = lib1990_obj@assays$RNA@counts, 
#			     path = "./data/lib1990_matrix",
#			     overwrite = TRUE)


# Demultiplexing and doublet detection
libs <- c("1984", "1988", "1990")

demuxlet_df <-
    "./demultiplexing/demuxlet/results/demuxlet_calls.tsv" |>
    read_tsv() |>
    rename("orig.ident" = "batch", "demuxlet_call" = "status") |>
    mutate(orig.ident = factor(orig.ident, levels = libs))

#scrublet_df <- 
#    glue::glue("./demultiplexing/scrublet/scrublet_calls_{libs}.tsv") |>
#    setNames(libs) |>
#    map_dfr(read_tsv, .id = "orig.ident") |>
#    rename("barcode" = "...1", 
#	   "scrublet_doublet_score" = "doublet_score") |>
#    mutate(orig.ident = factor(orig.ident, levels = libs)) |>
#    select(barcode, orig.ident, scrublet_doublet_score)

add_to_metadata <- function(x) {
    x@meta.data <- 
	x@meta.data |> 
	as_tibble(rownames = "barcode") |>
	left_join(demuxlet_df, join_by(barcode, orig.ident)) |>
	#left_join(scrublet_df, join_by(barcode, orig.ident)) |>
	column_to_rownames("barcode")
	
    x
} 

lib1984_obj <- add_to_metadata(lib1984_obj)
lib1988_obj <- add_to_metadata(lib1988_obj)
lib1990_obj <- add_to_metadata(lib1990_obj)

singlet_cells <- 
    list(lib1984_obj, lib1988_obj, lib1990_obj) |>
    map_dfr(function(x) x@meta.data |> 
	    as_tibble(rownames = "barcode")) |>
    filter(demuxlet_call == "SNG", 
	   #scrublet_doublet_score < 0.5,
	   dmm_type == "singlet") |>
    {function(x) split(x, x$orig.ident)}() |>
    map(~pull(., barcode))

lib1984_obj <- subset(lib1984_obj, cells = singlet_cells[["1984"]])
lib1988_obj <- subset(lib1988_obj, cells = singlet_cells[["1988"]])
lib1990_obj <- subset(lib1990_obj, cells = singlet_cells[["1990"]])

# Run PCA and Harmony to indentify outliers or contaminating cells
mt_ribo_genes <- 
    features |>
    filter(phenotype == "Gene Expression") |>
    filter(grepl("^MT-|^MRPS|^MRPL|^RPS|^RPL", gene_name))

bcells <- 
    merge(lib1984_obj, y = c(lib1988_obj, lib1990_obj), 
	  add.cell.ids = c("1984", "1988", "1990"))

bcells@meta.data$orig.ident <- paste0("BRI-", bcells@meta.data$orig.ident)

bcells <- bcells |>
    FindVariableFeatures() |>
    ScaleData(vars.to.regress = c("nCount_RNA", "percent_mt")) |>
    RunPCA()

set.seed(1L)
bcells <- bcells |>
    RunHarmony(group.by.vars = c("orig.ident", "donor_id"),
	       max_iter = 30,
	       reduction.save = "harmony")

sdev_plot <- 
    tibble(sdev = bcells@reductions$harmony@stdev) |>
    rownames_to_column("pc") |>
    mutate(pc = fct_inorder(pc)) |>
    ggplot(aes(pc, sdev)) +
    geom_line(aes(group = 1)) +
    scale_x_discrete(breaks = c(1, seq(5, 50, 5))) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    labs(x = "Principal component", y = "Standard deviation")

ggsave("./plots/v4_hpcs_sdev.png", sdev_plot, width = 4, height = 3)

# Umap and clustering
bcells <- bcells |>
    RunUMAP(reduction = "harmony", 
	    dims = 1:30,
	    seed.use = 1L,
	    reduction.name = "umap") |>
    FindNeighbors(dims = 1:30, reduction = "harmony", nn.eps = .5) |>
    FindClusters(resolution = 0.5)

Idents(bcells) <- "RNA_snn_res.0.5"

cluster_markers <- 
    FindAllMarkers(bcells, 
                   only.pos = TRUE,
                   min.pct = 0.1,
                   logfc.threshold = .5) |>
    as_tibble()

# Identify contaminating cell types and MALAT1-expressing cells
bad_clusters <-
    cluster_markers |>
    left_join(filter(features, phenotype == "Gene Expression"), 
	      join_by(gene == gene_id)) |>
    filter(gene_name %in% c("MALAT1", "GNLY", "S100A9")) |>
    group_by(gene_name) |>
    slice_max(avg_log2FC) |>
    pull(cluster)

bad_cells <- bcells@meta.data |>
    as_tibble(rownames = "barcode") |>
    filter(RNA_snn_res.0.5 %in% bad_clusters) |>
    pull(barcode)

bcells <- subset(bcells, cells = bad_cells, invert = TRUE)

bcells <- DietSeurat(bcells, dimreducs = NULL)
bcells@commands <- list()
bcells@meta.data$RNA_snn_res.0.5 <- NULL
bcells@meta.data$seurat_clusters <- NULL

Idents(bcells) <- "dmm_hto_call"

# Rerun analysis with the final object
# Scale and run PCAdd
bcells <- bcells |>
    FindVariableFeatures() |>
    ScaleData(vars.to.regress = c("nCount_RNA", "percent_mt")) |>
    RunPCA()

# Run Harmony correcting for batch
set.seed(1L)
bcells <- bcells |>
    RunHarmony(group.by.vars = c("orig.ident", "donor_id"),
	       max_iter = 30,
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
	   umap_1 = UMAP_1, umap_2 = UMAP_2) |>
    sample_frac(1L) |>
    mutate(barcode = fct_inorder(barcode),
	   hto = factor(hto, levels = names(stim_colors)),
	   cluster = factor(cluster, levels = sort(as.integer(levels(cluster)))))

write_tsv(umap_df, "./data/v4_umap_df.tsv")


umap_stim <-
    ggplot(umap_df, aes(umap_1, umap_2)) +
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

ggsave("plots/v4_umap_stim.png", umap_stim, width = 5, height = 4)

# Marker genes
Idents(bcells) <- "RNA_snn_res.0.5"

cluster_markers_2 <- 
    FindAllMarkers(bcells, 
                   only.pos = FALSE,
                   min.pct = 0.33,
                   logfc.threshold = .5) |>
    as_tibble() |>
    left_join(filter(features, phenotype == "Gene Expression") |> select(-phenotype), 
	      join_by(gene == gene_id))

write_tsv(cluster_markers_2, "./data/v4_cluster_markers.tsv")


# Save object
write_rds(bcells, "./data/seuratv4_qced.rds")


################################################################################
# Tests
################################################################################
stim_colors <- 
    read_tsv("../paper_plots/figure_colors.txt", col_names = c("stim", "time", "color")) |>
    unite("stim", c(stim, time), sep = " ") |>
    mutate(stim = paste0(stim, "h")) |>
    deframe()

bcells <- read_rds("./data/seuratv4_qced.rds")

adt_df <- 
    bcells@assays$ADT@data[c('IgD', 'CD27', 'CD69', 'CD86'), ] |>
    t() |>
    as_tibble(rownames = 'barcode')

umap_df <- 
    read_tsv("./data/v4_umap_df.tsv", col_types = "ccccfdd") |>
    mutate(hto = str_replace(hto, "IL4", "IL-4c"),
	   hto = str_replace(hto, "TLR7", "TLR7c"),
	   hto = str_replace(hto, "BCR", "BCRc"),
	   hto = str_replace(hto, "DN2", "DN2c")) |>
    left_join(adt_df, join_by(barcode)) |>
    mutate(hto = factor(hto, levels = names(stim_colors)))


igd_cd27_plot_density <- 
    ggplot(umap_df, aes(x = IgD, y = CD27)) +
    scale_fill_continuous(type = "viridis") +
    geom_vline(xintercept = .75, linetype = 2, color = "green") +
    geom_hline(yintercept = .25, linetype = 2, color = "green") +
    facet_wrap(~hto, nrow = 2) +
    geom_hex(aes(fill = after_stat(count/max(count))), bins = 100) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    guides(fill = "none")

umap_cd27 <-
    ggplot(data = umap_df |> arrange(CD27) |> mutate(barcode = fct_inorder(barcode)), 
	    aes(x = umap_1, y = umap_2)) +
    geom_point(aes(color = CD27), size = .5, stroke = 0) +
    scale_color_gradientn("CD27",
			  colors = c("grey85","#FFF7EC","#FEE8C8","#FDD49E","#FDBB84",
				     "#FC8D59","#EF6548","#D7301F","#B30000","#7F0000"),
			  label = function(x) sprintf("%.1f", x)) +
    theme_minimal() +
    theme(axis.title = element_blank(),
	  axis.text = element_blank(),
	  panel.grid = element_blank(),
	  ) +
    labs(color = NULL) +
    guides(color = guide_colorbar(barwidth = .5, barheight = 4))

umap_cd27_pos <- 
    ggplot(data = umap_df |> mutate(pos = CD27 >= 0.25) |> arrange(pos) |> mutate(barcode = fct_inorder(barcode)) ,
	    aes(x = umap_1, y = umap_2)) +
    geom_point(aes(color = pos), size = .5, stroke = 0) +
    scale_color_manual("CD27+", values = c("TRUE" = "midnightblue", "FALSE" = "grey85")) +
    theme_minimal() +
    theme(axis.title = element_blank(),
	  axis.text = element_blank(),
	  panel.grid = element_blank(),
	  ) +
    labs(color = NULL) +
    guides(color = guide_legend(override.aes = list(size = 3)))


set.seed(1)
umap_cd27_random <-
    ggplot(data = umap_df |> sample_frac(1) |> mutate(barcode = fct_inorder(barcode)), 
	    aes(x = umap_1, y = umap_2)) +
    geom_point(aes(color = CD27), size = .5, stroke = 0) +
    scale_color_gradientn("CD27",
			  colors = c("grey85","#FFF7EC","#FEE8C8","#FDD49E","#FDBB84",
				     "#FC8D59","#EF6548","#D7301F","#B30000","#7F0000"),
			  label = function(x) sprintf("%.1f", x)) +
    theme_minimal() +
    theme(axis.title = element_blank(),
	  axis.text = element_blank(),
	  panel.grid = element_blank(),
	  ) +
    labs(color = NULL) +
    guides(color = guide_colorbar(barwidth = .5, barheight = 4))

set.seed(1)
umap_cd27_pos_random <- 
    ggplot(data = umap_df |> mutate(pos = CD27 >= 0.25) |> sample_frac(1) |> mutate(barcode = fct_inorder(barcode)) ,
	    aes(x = umap_1, y = umap_2)) +
    geom_point(aes(color = pos), size = .5, stroke = 0) +
    scale_color_manual("CD27+", values = c("TRUE" = "midnightblue", "FALSE" = "grey85")) +
    theme_minimal() +
    theme(axis.title = element_blank(),
	  axis.text = element_blank(),
	  panel.grid = element_blank(),
	  ) +
    labs(color = NULL) +
    guides(color = guide_legend(override.aes = list(size = 3)))

#library(patchwork)


subsets_df <- 
    adt_df |>
    mutate(
	   "IgD<sup>+</sup>  CD27<sup>+</sup>"  = ifelse(IgD >= .75 & CD27 >= .25, 1, 0),
	   "IgD<sup>–</sup>  CD27<sup>–</sup>" = ifelse(IgD < .75 & CD27 < .25, 1, 0),
	   "IgD<sup>–</sup>  CD27<sup>+</sup>" = ifelse(IgD < .75 & CD27 >= .25, 1, 0),
	   "IgD<sup>+</sup>  CD27<sup>–</sup>" = ifelse(IgD >= .75 & CD27 < .25, 1, 0)
	   ) |>
    select(barcode, contains("sup")) |>
    pivot_longer(-barcode, names_to = 'subtype') |>
    mutate(subtype = fct_inorder(subtype)) |>
    left_join(select(umap_df, barcode, hto)) |>
    group_by(hto, subtype) |>
    summarise(pct = mean(value)) |>
    ungroup()

proportions_plot <- 
    ggplot(subsets_df, aes(x = hto, y = pct)) +
    geom_col(aes(fill = hto)) +
    scale_fill_manual(values = stim_colors) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, .25, .5, .75, 1)) +
    facet_wrap(~subtype, nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
	  panel.grid.major.x = element_blank(),
	  panel.grid.minor.y = element_blank(),
	  panel.spacing.x = unit(1, "lines"),
	  axis.ticks.x = element_blank(),
	  strip.text = ggtext::element_markdown(size = 8),
	  strip.clip = "off",
	  plot.background = element_rect(color = "white", fill = "white")
	  ) +
    guides(fill = "none") +
    labs(x = NULL, y = NULL)


library(cowplot)


fig_a_title <- 
    ggdraw() + 
    draw_label(
	       "Gating",
	       x = 0,
	       size = 9,
	       hjust = 0
	       ) +
    theme(text = element_text(size = 12),
	  plot.margin = margin(l = 1.25, unit = "lines"))

fig_a <- 
    plot_grid(fig_a_title,
	      igd_cd27_plot_density, 
	      ncol = 1, 
	      rel_heights = c(.1, 1),
	      labels = 'a', label_size = 12, vjust = .75)

fig_b_title <- 
    ggdraw() + 
    draw_label(
	       "UMAP with droplets sorted by CD27 expression",
	       x = 0,
	       size = 9,
	       hjust = 0
	       ) +
    theme(text = element_text(size = 12),
	  plot.margin = margin(l = 1.25, unit = "lines"))

fig_b <- 
    plot_grid(fig_b_title,
	      plot_grid(umap_cd27, umap_cd27_pos),
	      ncol = 1, 
	      rel_heights = c(.1, 1),
	      labels = 'b', label_size = 12, vjust = .75)

fig_c_title <- 
    ggdraw() + 
    draw_label(
	       "UMAP with droplets randomly sorted",
	       x = 0,
	       size = 9,
	       hjust = 0
	       ) +
    theme(text = element_text(size = 12),
	  plot.margin = margin(l = 1.25, unit = "lines"))

fig_c <- 
    plot_grid(fig_c_title,
	      plot_grid(umap_cd27_random, umap_cd27_pos_random),
	      ncol = 1, 
	      rel_heights = c(.1, 1),
	      labels = 'c', label_size = 12, vjust = .75)

fig_d_title <- 
    ggdraw() + 
    draw_label(
	       "Proportion of cell subsets in each condition",
	       x = 0,
	       size = 9,
	       hjust = 0
	       ) +
    theme(text = element_text(size = 12),
	  plot.margin = margin(l = 1.25, unit = "lines"))

fig_d <- 
    plot_grid(fig_d_title,
	      proportions_plot,
	      ncol = 1, 
	      rel_heights = c(.1, 1),
	      labels = 'd', label_size = 12, vjust = .75)

plot_out <- 
    plot_grid(fig_a, fig_b, fig_c, fig_d, ncol = 1, rel_heights = c(1, .8, .8, .5)) +
    theme(plot.background = element_rect(color = "white", fill = "white"))


ggsave("./plots/bcell_cd27.png", plot_out, width = 6.5, height = 8.5)

# Try kmeans
cd27_vec <- 
    adt_df |>
    select(barcode, CD27) |>
    deframe()

head(cd27_vec)
km <- kmeans(cd27_vec, centers = 2)

km_df <- km$cluster |> enframe("barcode", "k")

umap_df_k <- 
    left_join(umap_df, km_df, join_by(barcode))

set.seed(1)
umap_cd27_kmeans <- 
    ggplot(data = umap_df_k |> sample_frac(1) |> mutate(barcode = fct_inorder(barcode)),
	    aes(x = umap_1, y = umap_2)) +
    geom_point(aes(color = factor(k)), size = .5, stroke = 0) +
    scale_color_manual("CD27+", values = c("2" = "midnightblue", "1" = "grey85")) +
    theme_bw() +
    theme(axis.title = element_blank(),
	  axis.text = element_blank(),
	  panel.grid = element_blank(),
	  ) +
    labs(color = NULL) +
    guides(color = guide_legend(override.aes = list(size = 3)))


kmean_violin <- 
    ggplot(umap_df_k, aes(x = factor(k), y = CD27)) +
    geom_violin() +
    theme_bw()

ggsave("./plots/umap_cd27_kmeans.png", 
       plot_grid(kmean_violin, umap_cd27_kmeans, nrow = 1, align = "h", rel_widths = c(.6, 1)), 
       width = 8, height = 4)


library(mclust)

fit <- Mclust(cd27_vec[cd27_vec>0], G = 2)

# Sample points between the peaks
low_peak <- min(fit$parameters$mean)
high_peak <- max(fit$parameters$mean)
test_points <- seq(low_peak, high_peak, length.out = 100)

# Find where density is lowest
densities <- dens(modelName = fit$modelName, 
                  data = test_points,
                  parameters = fit$parameters)

threshold <- test_points[which.min(densities)]

cd27_hist <- 
    ggplot(umap_df |> filter(CD27 > 0), aes(x = CD27)) +
    geom_histogram(bins = 100) +
    geom_vline(xintercept = c(low_peak, high_peak)) +
    geom_vline(xintercept = threshold, color = "red", linetype = 3) +
    scale_x_continuous() +
    theme_minimal() +
    theme(plot.background = element_rect(color = "white", fill = "white"))

ggsave("./plots/cd27_hist.png", cd27_hist, width = 4, height = 1.5)
