# Set large RAM limit for R
unix::rlimit_as(1e12)

library(Seurat)
library(SeuratDisk)
library(harmony)
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

# Import Seurat objects
lib1984 <- readRDS("./data/seurat_1984_qced.rds")
lib1988 <- readRDS("./data/seurat_1988_qced.rds")
lib1990 <- readRDS("./data/seurat_1990_qced.rds")

DefaultAssay(lib1984) <- DefaultAssay(lib1988) <- DefaultAssay(lib1990) <- "RNA"

# Merge datasets
bcells <- 
    merge(lib1984, y = c(lib1988, lib1990), add.cell.ids = c("1984", "1988", "1990"))

# Rename genes
## ensure that genes are in the same order
all(rownames(bcells[["RNA"]]@data) == genes_df$gene_id)

## Rename
rownames(bcells[["RNA"]]@counts) <- genes_df$gene_lab
rownames(bcells[["RNA"]]@data) <- genes_df$gene_lab
rownames(bcells[["RNA"]]@meta.features) <- genes_df$gene_lab

# Scale and run PCA
bcells <- bcells |>
    {function(x) ScaleData(x, features = rownames(x))}() |>
    FindVariableFeatures() |>
    {function(x) RunPCA(x, features = VariableFeatures(x))}()

# Run Harmony
# Correct for batch
set.seed(1L)
bcells <- bcells |>
    RunHarmony(group.by.vars = "orig.ident", 
	       max.iter.harmony = 30,
	       reduction.save = "harmony")

# UMAP
bcells <- bcells |>
    RunUMAP(reduction = "harmony", 
	    dims = 1:35,
	    seed.use = 1L,
	    reduction.name = "umap")

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

# Reduce Seurat object by keeping only essential elements for visualization
bcells_out <- 
    DietSeurat(bcells,
	       assays = "RNA",
	       dimreducs = "umap",
	       counts = FALSE,
	       scale.data = FALSE)

# Save data
SaveH5Seurat(bcells_out, "bcells_harmony.h5Seurat", overwrite = TRUE)
Convert("bcells_harmony.h5Seurat", dest = "h5ad", overwrite = TRUE)

