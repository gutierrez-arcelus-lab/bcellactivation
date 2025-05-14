library(tidyverse)
library(vroom)
library(Seurat)
library(hdf5r)


# Bulk RNA-seq
dat <- read_rds("./deseq_normalized_counts.rds")

write_lines(unique(dat$gene_label), "./bulk_genes.txt")


# Single-cell data
sc_data <- SeuratDisk::LoadH5Seurat("../../citeseq/data/bcells.h5Seurat")

genes_expressed <- 
    GetAssayData(object = sc_data, assay = "RNA", slot = "data") |>
    apply(1, function(x) sum(x > 0)) |>
    {function(x) which(x >= 10)}() |>
    names()

sc_data_sub <- subset(sc_data, features = genes_expressed)

expr_mat <- 
    as.matrix(GetAssayData(sc_data_sub, assay = "RNA", slot = "data"))

chunk_dims <- c(1, ncol(expr_mat)) 

h5file <- H5File$new("./singlecell/bcells_expressed.h5", mode = "w")

dset <- 
    h5file$create_dataset(
			  name = "expr_data",
			  dtype = h5types$H5T_NATIVE_FLOAT,  
			  dims = dim(expr_mat),
			  chunk_dims = chunk_dims
)

for (i in 1:nrow(expr_mat)) 
  dset[i, ] <- expr_mat[i, , drop = FALSE]

h5file$close_all()

write_lines(rownames(expr_mat), "./singlecell/genes.txt")
write_lines(colnames(expr_mat), "./singlecell/cells.txt")

sc_data_sub@meta.data |>
    as_tibble(rownames = "barcode") |>
    write_tsv("./singlecell/metadata.tsv")

Embeddings(sc_data_sub, "umap") |>
    as_tibble(rownames = "barcode") |>
    write_tsv("./singlecell/umap_data.tsv")

