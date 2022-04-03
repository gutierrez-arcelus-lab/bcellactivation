library(tidyverse)
library(rhdf5)

# Examine scDRS example h5 data
h5ls("./data/expr.h5ad")

xdata <- h5read("./data/expr.h5ad", "/X")
class(xdata)                            
dim(xdata)

xdata[1:5, 1:10]

cellid <- h5read("./data/expr.h5ad", "/obs/cell_id")
class(cellid)
dim(cellid)
head(cellid)

genes <- h5read("./data/expr.h5ad", "/var/gene") 
head(genes)

# B cells
bcells_singlet <- read_rds("./data/bcells_singlet_seurat.rds")

counts <- as.matrix(bcells_singlet@assays$RNA@counts)
bcell_ids <- array(colnames(counts))
bcell_genes <- as.array(rownames(counts))
dimnames(counts) <- NULL

h5out <- "./data/bcells_singlet_seurat.h5ad"
unlink(h5out)
h5createFile(h5out)

h5createDataset(h5out, "X", 
		dims = dim(counts), 
		storage.mode = "double")
		
h5write(counts, h5out, "X")

h5createGroup(h5out, "obs")
h5write(bcell_ids, h5out, name = "/obs/cell_id")

h5createGroup(h5out, "var")
h5write(bcell_genes, h5out, name = "/var/gene")

h5closeAll()

h5ls(h5out)
#h5read(h5out, "/var/gene")
