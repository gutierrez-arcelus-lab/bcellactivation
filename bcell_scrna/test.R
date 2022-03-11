library(tidyverse)
library(Seurat)
library(ggbeeswarm)


cellranger_dir <- 
    file.path("/lab-share/IM-Gutierrez-e2/Public/scRNA/SN0231064/KW9100_Maria",
              "210726_10X_KW9100-2_bcl/cellranger-6.0.1/GRCh38/BRI-1283/outs",
              "filtered_feature_bc_matrix")

features_df <- file.path(cellranger_dir, "features.tsv.gz") %>%
    read_tsv(col_names = c("gene_id", "gene_name", "phenotype"))

mt_genes <- features_df %>%
    filter(phenotype == "Gene Expression", 
           grepl("^MT-", gene_name)) %>%
    pull(gene_id)

ribo_genes <- features_df %>%
    filter(phenotype == "Gene Expression", 
           grepl("^RPS\\d+|^RPL\\d+", gene_name))

data10x <- Read10X(cellranger_dir, gene.column = 1)

gene_exp <- data10x[["Gene Expression"]]

antibody <- data10x[["Antibody Capture"]] %>%
    .[!grepl("^Hashtag", rownames(.)), ]

hashtags <- data10x[["Antibody Capture"]] %>%
    .[grepl("^Hashtag", rownames(.)), ]


bcells <- CreateSeuratObject(counts = gene_exp, project = "bcells",
                             min.cells = 3, min.features = 200)

bcells[["percent.mt"]] <- PercentageFeatureSet(bcells, features = mt_genes)

bcells_meta <- bcells@meta.data %>%
    rownames_to_column("barcode") %>%
    pivot_longer(-(barcode:orig.ident), names_to = "feature")

ggplot(bcells_meta, aes(orig.ident, value)) +
    geom_jitter(alpha = .25, size = .25) +
    geom_violin(alpha = .25) +
    facet_wrap(~feature, nrow = 1, scales = "free") +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank())

n_genes <- Matrix::colSums(gene_exp > 0)

bcells_mt_ngenes <- bcells_meta %>%
    filter(feature == "percent.mt") %>%
    pivot_wider(names_from = feature, values_from = value) %>%
    left_join(tibble(barcode = names(n_genes),
                     n_genes = n_genes), 
              by = "barcode")

ggplot(bcells_mt_ngenes, aes(percent.mt, n_genes)) +
    geom_jitter(size = .25, alpha = .25) +
    geom_density2d_filled(alpha = .5) +
    geom_density_2d(size = .25, color = "black") +
    geom_vline(xintercept = c(5, 10, 20), linetype = 2) +
    scale_x_continuous(breaks = c(0, 5, 10, 20, 40, 60, 80)) +
    scale_y_continuous(breaks = seq(0, 1e4, by = 1e3)) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          legend.position = "none") +
    labs(x = "% Mitochondrial reads", y = "# of genes")



bcells <- bcells %>%
    subset(subset = nFeature_RNA > 500 & percent.mt < 20) %>%
    NormalizeData()
    
bcells <- FindVariableFeatures(bcells,
                               selection.method = "vst",
                               nfeatures = 2000)
    
bcells <- ScaleData(bcells, features = rownames(bcells))

bcells <- RunPCA(bcells, features = VariableFeatures(bcells))

ElbowPlot(bcells)

bcells <- FindNeighbors(bcells, dims = 1:20)
bcells <- FindClusters(bcells, resolution = .2)
bcells <- RunUMAP(bcells, dims = 1:20)

DimPlot(bcells, reduction = "umap")

as.data.frame(bcells@reductions$umap@cell.embeddings) %>%
    as_tibble() %>%
    mutate(cluster = bcells@meta.data$seurat_clusters) %>%
    ggplot(aes(UMAP_1, UMAP_2, color = cluster)) +
    geom_point() +
    scale_color_viridis_d() +
    theme_bw()






    
    