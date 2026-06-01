# ==============================================================================
# Description:  Core processing script for the CITE-seq data. Handles multi-modal 
#               data import (RNA/ADT/HTO), quality control, HTO-based demultiplexing 
#               (demuxmix), integration with SNP-based demultiplexing (demuxlet), 
#               contaminant removal, Harmony batch correction, UMAP projection, 
#               and cluster marker identification.
# Input:        1. ../data/cellranger/.../filtered_feature_bc_matrix (Cell Ranger)
#               2. ../1-demultiplexing/demuxlet/results/demuxlet_calls.tsv
#               3. ../../paper_figures/figure_colors_final.txt (Plot styling)
# Output:       1. ./data/v4_seurat_qced.rds (Final, clean Seurat object)
#               2. ./data/v4_umap_df.tsv (UMAP coordinates)
#               3. ./data/v4_cluster_markers.tsv (Cluster-specific genes)
#               4. ./plots/v4_hpcs_sdev.png & v4_umap_stim.png
# ==============================================================================

# ==============================================================================
# 1. Environment Setup & Helper Functions
# ==============================================================================

# Increase memory limit for handling large single-cell matrices
unix::rlimit_as(1e12)

library(Seurat)
library(demuxmix)
library(harmony)
library(tidyverse)

# Ensure output directories exist
if (!file.exists("data")) dir.create("data")
if (!file.exists("plots")) dir.create("plots")

# Import standardized publication colors for the stimulation conditions
stim_colors <- 
    "../../paper_figures/figure_colors_final.txt" |>
    read_tsv() |>
    mutate(stim = glue::glue("{Condition} {Time}h")) |>
    select(stim, Hex) |>
    deframe()

# Function to import 10x multi-modal data and initialize Seurat object
make_seurat <- function(cellranger_path, project_id, hto_names = NULL, mito_ids, ribo_ids) {

    # Read the 10x hdf5/matrix files (using gene IDs as row names)
    data10x <- Read10X(cellranger_path, gene.column = 1)

    # Separate ADT (surface proteins) from HTO (hashing antibodies)
    antibody_mtx <- data10x[["Antibody Capture"]] |>
        {function(x) x[!grepl("^Hashtag", rownames(x)), ]}()

    if (!is.null(hto_names)) { 

	hashtags_mtx <- data10x[["Antibody Capture"]] |>
	    {function(x) x[grepl("^Hashtag", rownames(x)), ]}()

	rownames(hashtags_mtx) <- setNames(hto_names[rownames(hashtags_mtx)], NULL)
    }

    # Initialize Seurat with RNA assay and LogNormalize
    seuratobj <- 
	CreateSeuratObject(counts = data10x[["Gene Expression"]], project = project_id) |>
        NormalizeData(assay = "RNA", normalization.method = "LogNormalize")

    # Add ADT assay and apply Centered Log Ratio (CLR) normalization
    if ( nrow(antibody_mtx) > 0 ) {
	
	rownames(antibody_mtx) <- sub("_prot$", "", rownames(antibody_mtx))
	
	seuratobj[["ADT"]] <- CreateAssayObject(counts = antibody_mtx)
    
	seuratobj <- NormalizeData(seuratobj, assay = "ADT", normalization.method = "CLR", margin = 2)
    }
   
    # Add HTO assay and apply CLR normalization
    if (!is.null(hto_names)) {
	
	seuratobj[["HTO"]] <- CreateAssayObject(counts = hashtags_mtx)

	seuratobj <- NormalizeData(seuratobj, assay = "HTO", normalization.method = "CLR", margin = 2)
    }

    # Calculate mitochondrial and ribosomal percentages for QC
    seuratobj[["percent_mt"]] <- PercentageFeatureSet(seuratobj, features = mito_ids)
    seuratobj[["percent_ribo"]] <- PercentageFeatureSet(seuratobj, features = ribo_ids)

    seuratobj
}

# Function to demultiplex hashing antibodies using demuxmix
run_demuxmix <- function(x) {
   
    hto <- as.matrix(x@assays$HTO@counts)

    # Provide RNA feature counts to demuxmix to improve probability estimates
    rna <- x@meta.data |> 
	select(nFeature_RNA) |>
	as_tibble(rownames = "barcode") |>
	deframe()

    rna <- rna[colnames(hto)]
    
    # Run the demuxmix negative binomial mixture model
    dmm <- demuxmix(hto, rna = rna, maxIter = 500, tol = 1e-4)
    
    classes <- 
	dmmClassify(dmm) |>
	as_tibble(rownames = "barcode") |>
	select(barcode, dmm_hto_call = HTO, dmm_prob = Prob, dmm_type = Type)
    
    # Append demultiplexing results to Seurat metadata
    x@meta.data <- 
	x@meta.data |>
	as_tibble(rownames = "barcode") |>
	left_join(classes, join_by(barcode)) |>
	column_to_rownames("barcode")

    x
}

# Function to append the SNP-based demuxlet calls to Seurat metadata
add_to_metadata <- function(x) {
    x@meta.data <- 
	x@meta.data |> 
	as_tibble(rownames = "barcode") |>
	left_join(demuxlet_df, join_by(barcode, orig.ident)) |>
	#left_join(scrublet_df, join_by(barcode, orig.ident)) |>
	column_to_rownames("barcode")
	
    x
} 

# ==============================================================================
# 2. Data Initialization
# ==============================================================================
lib1984_dir <- "../data/cellranger/1984/filtered_feature_bc_matrix"
lib1988_dir <- "../data/cellranger/1988/filtered_feature_bc_matrix"
lib1990_dir <- "../data/cellranger/1990/filtered_feature_bc_matrix"

# Feature IDs
features <- 
    file.path(lib1984_dir, "features.tsv.gz") |>
    read_tsv(col_names = c("gene_id", "gene_name", "phenotype"))

# Identify mitochondrial and ribosomal genes across the dataset
mt_genes <- features |>
    filter(grepl("^MT-", gene_name)) |>
    pull(gene_id)

ribo_genes <- features |>
    filter(grepl("^RPS\\d+|^RPL\\d+", gene_name)) |>
    pull(gene_id)

# Map hashtag oligonucleotide (HTO) IDs to the physical stimulation conditions
hashtags <- 
    c("Hashtag6" = "Unstim 0h",
      "Hashtag7" = "IL-4c 24h",
      "Hashtag8" = "IL-4c 72h",
      "Hashtag9" = "BCRc 24h",
      "Hashtag10" = "BCRc 72h",
      "Hashtag12" = "TLR7c 24h",
      "Hashtag13" = "TLR7c 72h",
      "Hashtag14" = "DN2c 24h",
      "Hashtag15" = "DN2c 72h")

# Create Seurat objects and run HTO demultiplexing
lib1984_obj <- 
    make_seurat(lib1984_dir, 
		project_id = "1984", 
		hto_names = hashtags,
		mito_ids = mt_genes, ribo_ids = ribo_genes) |>
    run_demuxmix()

lib1988_obj <- 
    make_seurat(lib1988_dir, 
		project_id = "1988", 
		hto_names = hashtags,
		mito_ids = mt_genes, ribo_ids = ribo_genes) |>
    run_demuxmix()

lib1990_obj <- 
    make_seurat(lib1990_dir, 
		project_id = "1990", 
		hto_names = hashtags,
		mito_ids = mt_genes, ribo_ids = ribo_genes) |>
    run_demuxmix()

# ==============================================================================
# 3. Quality Control & Contaminant Filtering
# ==============================================================================

# A. Basic QC filtering: Keep cells with > 500 genes and < 10% mitochondrial reads
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

# B. Integration of SNP Demultiplexing (demuxlet)
libs <- c("1984", "1988", "1990")

demuxlet_df <-
    "../1-demultiplexing/demuxlet/results/demuxlet_calls.tsv" |>
    read_tsv() |>
    rename("orig.ident" = "batch", "demuxlet_call" = "status") |>
    mutate(orig.ident = factor(orig.ident, levels = libs)) |>
    add_count(orig.ident, donor_id) |>
    mutate(demuxlet_call = case_when(demuxlet_call == "SNG" & n < 50 ~ "ERROR",
				     .default = demuxlet_call)) |>
    select(-n)

lib1984_obj <- add_to_metadata(lib1984_obj)
lib1988_obj <- add_to_metadata(lib1988_obj)
lib1990_obj <- add_to_metadata(lib1990_obj)

# Keep confirmed singlets (Singlet by both genetics and hashtagging)
singlet_cells <- 
    list("1984" = lib1984_obj, "1988" = lib1988_obj, "1990" = lib1990_obj) |>
    map(function(x) x@meta.data |> 
	as_tibble(rownames = "barcode") |>
	filter(demuxlet_call == "SNG", dmm_type == "singlet") |>
	pull(barcode)
    )

lib1984_obj <- subset(lib1984_obj, cells = singlet_cells[["1984"]])
lib1988_obj <- subset(lib1988_obj, cells = singlet_cells[["1988"]])
lib1990_obj <- subset(lib1990_obj, cells = singlet_cells[["1990"]])

# ==============================================================================
# 4. Pre-processing & Removal of Ambient Contaminants
# ==============================================================================

# Merge the three batches into a single object
bcells <- 
    merge(lib1984_obj, y = c(lib1988_obj, lib1990_obj), 
	  add.cell.ids = c("1984", "1988", "1990"))

bcells@meta.data$orig.ident <- paste0("BRI-", bcells@meta.data$orig.ident)

# Initial round of scaling and dimensionality reduction
bcells <- bcells |>
    FindVariableFeatures() |>
    ScaleData(vars.to.regress = c("nCount_RNA", "percent_mt")) |>
    RunPCA()

# Run Harmony. Grouping by both orig.ident (library batch) and donor_id accounts 
# for both technical sequencing variation and baseline genetic variation.
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

# Cluster the data
bcells <- bcells |>
    RunUMAP(reduction = "harmony", 
	    dims = 1:30,
	    seed.use = 1L,
	    reduction.name = "umap") |>
    FindNeighbors(dims = 1:30, reduction = "harmony", nn.eps = .5) |>
    FindClusters(resolution = 0.5)

Idents(bcells) <- "RNA_snn_res.0.5"

# Find broad cluster markers
cluster_markers <- 
    FindAllMarkers(bcells, 
                   only.pos = TRUE,
                   min.pct = 0.1,
                   logfc.threshold = .5) |>
    as_tibble()

# Identify and remove clusters defined by 
# MALAT1 (dying cells), 
# GNLY (NK/T cells), 
# or S100A9 (Monocytes/Neutrophils).
bad_clusters <-
    cluster_markers |>
    left_join(filter(features, phenotype == "Gene Expression"), 
	      join_by(gene == gene_id)) |>
    filter(gene_name %in% c("MALAT1", "GNLY", "S100A9")) |>
    group_by(gene_name) |>
    slice_max(avg_log2FC) |>
    pull(cluster)

bad_cells <- 
    bcells@meta.data |>
    as_tibble(rownames = "barcode") |>
    filter(RNA_snn_res.0.5 %in% bad_clusters) |>
    pull(barcode)

bcells <- subset(bcells, cells = bad_cells, invert = TRUE)

# Strip out old reductions and metadata columns before the final processing run
bcells <- DietSeurat(bcells, dimreducs = NULL)
bcells@commands <- list()
bcells@meta.data$RNA_snn_res.0.5 <- NULL
bcells@meta.data$seurat_clusters <- NULL

Idents(bcells) <- "dmm_hto_call"

# ==============================================================================
# 5. Final Processing of Cleaned B Cells
# ==============================================================================

# Re-run the entire scaling, PCA, Harmony, and clustering pipeline on the clean cells
bcells <- 
    bcells |>
    FindVariableFeatures() |>
    ScaleData(vars.to.regress = c("nCount_RNA", "percent_mt")) |>
    RunPCA()

# Run Harmony correcting for batch
set.seed(1L)
bcells <- 
    bcells |>
    RunHarmony(group.by.vars = c("orig.ident", "donor_id"),
	       max_iter = 30,
	       reduction.save = "harmony") |>
    FindNeighbors(dims = 1:30, reduction = "harmony") |>
    FindClusters(resolution = 0.5) |>
    RunUMAP(reduction = "harmony", 
	    dims = 1:30,
	    seed.use = 1L,
	    reduction.name = "umap")

# Generate UMAP visualization metadata
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

# ==============================================================================
# 6. Final Cluster Marker Identification
# ==============================================================================
Idents(bcells) <- "RNA_snn_res.0.5"

cluster_markers_2 <- 
    FindAllMarkers(bcells, 
                   only.pos = FALSE,
                   min.pct = 0.33,
                   logfc.threshold = .5) |>
    as_tibble() |>
    left_join(filter(features, phenotype == "Gene Expression") |> select(-phenotype), 
	      join_by(gene == gene_id))

# ==============================================================================
# 7. Data Export
# ==============================================================================
write_rds(bcells, "./data/v4_seurat_qced.rds")
write_tsv(umap_df, "./data/v4_umap_df.tsv")
write_tsv(cluster_markers_2, "./data/v4_cluster_markers.tsv")
