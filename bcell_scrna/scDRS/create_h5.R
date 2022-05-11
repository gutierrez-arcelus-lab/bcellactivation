library(tidyverse)
library(reticulate)
library(Seurat)
library(SingleCellExperiment)

seurat <- read_rds("./data/expression/bcells_singlet_seurat.rds")

sc <- import("scanpy")

anndata <- sc$AnnData(
    X = as.matrix(t(GetAssayData(seurat, slot = "data"))),
    obs = seurat[[]],
    var = GetAssay(seurat)[[]]
)

sc$AnnData$write_h5ad(anndata, "./data/expression/bcells_singlet_seurat.h5ad")
