library(tidyverse)
library(reticulate)
library(Seurat)
library(SingleCellExperiment)

seurat <- read_rds("./data/bcells_singlet_seurat.rds")

str(seurat)

sc <- import("scanpy")

anndata <- sc$AnnData(
    X = as.matrix(t(GetAssayData(seurat, slot = "data"))),
    obs = seurat[[]],
    var = GetAssay(seurat)[[]]
)

anndata

sc$AnnData$write_h5ad(anndata, "./data/bcells_singlet_seurat.h5ad")


