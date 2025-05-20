# Set large RAM limit for R
unix::rlimit_as(1e12)

library(Seurat)
library(SeuratDisk)
library(dplyr)
library(readr)

# Gene ids
cellranger_dir_1984 <- 
    file.path("/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/B_cells_citeseq",
	      "SN0268787/broad/hptmp/curtism/bwh10x/KW10456_Maria",
	      "221101_10X_KW10456_bcl/cellranger-6.1.1/GRCh38/BRI-1984/outs",
	      "filtered_feature_bc_matrix")

features_df <- 
    file.path(cellranger_dir_1984, "features.tsv.gz") |>
    read_tsv(col_names = c("gene_id", "gene_name", "phenotype"))

# append gene IDs to disambiguate genes with the same gene name
genes_df <- 
    features_df |> 
    as_tibble() |>
    filter(phenotype == "Gene Expression") |>
    select(gene_id, gene_name) |>
    add_count(gene_name) |>
    mutate(gene_lab = ifelse(n == 1, gene_name, paste(gene_name, gene_id, sep = "-"))) |>
    select(gene_id, gene_lab)

# Import Seurat object
bcells <- read_rds("./data/seuratv4_qced.rds") 

# Rename genes
## ensure that genes are in the same order
all(rownames(bcells[["RNA"]]@data) == genes_df$gene_id)

## then rename
rownames(bcells[["RNA"]]@counts) <- genes_df$gene_lab
rownames(bcells[["RNA"]]@data) <- genes_df$gene_lab
rownames(bcells[["RNA"]]@meta.features) <- genes_df$gene_lab

# Append ADT data to RNA
## First add the 'prot' suffix to ADT protein names
rownames(bcells[['ADT']]@data) <- 
    paste0("prot_", rownames(bcells[['ADT']]@data))

# ensure that cells are in the same order
all(colnames(bcells[['RNA']]@data) == colnames(bcells[['ADT']]))

# append data
bcells[['RNA']]@data <- rbind(bcells[['RNA']]@data, bcells[['ADT']]@data)

adt_meta_features <- 
    matrix(NA, 
	   nrow = nrow(bcells[['ADT']]@data),
	   ncol = ncol(bcells[['RNA']]@meta.features)) |>
    as.data.frame() |>
    setNames(colnames(bcells[['RNA']]@meta.features)) |>
    `rownames<-`(rownames(bcells[['ADT']]@data))

bcells[['RNA']]@meta.features <- 
    rbind(bcells[['RNA']]@meta.features, adt_meta_features)

# update var features
bcells$RNA@var.features <- 
    genes_df |> 
    filter(gene_id %in% bcells$RNA@var.features) |>
    mutate(gene_id = factor(gene_id, levels = bcells$RNA@var.features)) |>
    arrange(gene_id) |>
    pull(gene_lab)

# Reduce Seurat object by keeping only essential elements for visualization
bcells_out <- 
    DietSeurat(bcells,
	       assays = "RNA",
	       dimreducs = "umap",
	       counts = FALSE,
	       graphs = FALSE,
	       scale.data = FALSE,
	       misc = FALSE)

bcells_out@commands <- list()

bcells_out@meta.data <- 
    bcells_out@meta.data |>
    select(-RNA_snn_res.0.5) |>
    mutate(seurat_clusters = paste0("C", seurat_clusters))

# Save data
SaveH5Seurat(bcells_out, "./data/bcells.h5Seurat", overwrite = TRUE)
Convert("./data/bcells.h5Seurat", dest = "h5ad", overwrite = TRUE)
