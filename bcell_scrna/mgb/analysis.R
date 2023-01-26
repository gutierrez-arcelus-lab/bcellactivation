# single-cell data analysis
library(Seurat)
library(miQC)
library(scater)

# Data wrangling
library(tidyverse)

# Plotting
library(tidytext)
library(ggridges)
library(RColorBrewer)
library(scico)
library(ggsci)
library(cowplot)


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

cellranger_dir <- 
    file.path("/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/B_cells_citeseq",
	      "SN0268787/broad/hptmp/curtism/bwh10x/KW10456_Maria",
	      "221101_10X_KW10456_bcl/cellranger-6.1.1/GRCh38/BRI-1984/outs",
	      "filtered_feature_bc_matrix")

features_df <- file.path(cellranger_dir, "features.tsv.gz") |>
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

data10x <- Read10X(cellranger_dir, gene.column = 1)

antibody_mtx <- data10x[["Antibody Capture"]] |>
    {function(x) x[!grepl("^Hashtag", rownames(x)), ]}()

rownames(antibody_mtx) <- rownames(antibody_mtx) |>
    sub("_prot$", "", x = _) |>
    gsub("_", ".", x = _)

hashtags_mtx <- data10x[["Antibody Capture"]] |>
    {function(x) x[grepl("^Hashtag", rownames(x)), ]}()

rownames(hashtags_mtx) <- setNames(stims[rownames(hashtags_mtx)], NULL)

stim_colors <- c("grey80", "grey40", "black", brewer.pal(n = 6, "Paired"))
names(stim_colors) <- stims

# Create object
bcells <- CreateSeuratObject(counts = data10x[["Gene Expression"]], project = "MGB")
bcells[["ADT"]] <- CreateAssayObject(counts = antibody_mtx)
bcells[["HTO"]] <- CreateAssayObject(counts = hashtags_mtx)

bcells <- bcells |>
    NormalizeData(normalization.method = "LogNormalize", margin = 1) |>
    NormalizeData(assay = "HTO", normalization.method = "CLR", margin = 1) |>
    NormalizeData(assay = "ADT", normalization.method = "CLR", margin = 2)

bcells[["percent_mt"]] <- PercentageFeatureSet(bcells, features = mt_genes)
bcells[["percent_ribo"]] <- PercentageFeatureSet(bcells, features = ribo_genes)

# Demuxlet
demuxlet_df <- read_tsv("./demuxlet/demuxlet_results.best") |>
    select(barcode = BARCODE, best = BEST) |>
    extract(best, c("status", "sample"), "([^-]+)-(.+)")

demuxlet_df |>
    filter(status == "SNG") |>
    count(sample)


# HTO counts
hto <- as_tibble(t(bcells@assays$HTO@counts), rownames = "barcode") |>
    pivot_longer(-barcode, names_to = "stim") |>
    group_by(barcode) |>
    mutate(hto_max = stim[which.max(value)]) |>
    ungroup() |>
    mutate_at(vars(stim, hto_max), ~factor(., levels = stims))

png("./plots/test.png", units = "in", res = 300, width = 6, height = 3)
ggplot(hto, aes(x = log10(value + 1))) +
    geom_density(aes(fill = stim), size = .25, alpha = .5) +
    scale_fill_manual(values = stim_colors) +
    facet_wrap(~hto_max, nrow = 2) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(x = "log10 (counts + 1)")
dev.off()

# QC
bcells_meta <- bcells@meta.data |>
    as_tibble(rownames = "barcode") 

qc_plot <- bcells_meta |>
    select(barcode, n_genes = nFeature_RNA, percent_mt) |> 
    ggplot(aes(x = n_genes, y = percent_mt)) +
    geom_point(size = .25, alpha = .2) +
    geom_vline(xintercept = 700, linetype = 2, color = "tomato3") +
    geom_hline(yintercept = 7, linetype = 2, color = "tomato3") +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    labs(x = "Number of genes", y = "% reads from Mitochondria")

ggsave("./plots/test.png", qc_plot, width = 4, height = 3)
   

demuxlet_singlets <- demuxlet_df |>
    filter(status == "SNG")

cells_keep <- bcells_meta |>
    filter(nFeature_RNA > 700, percent_mt <= 7) |>
    filter(barcode %in% demuxlet_singlets$barcode)

bcells_pass <- subset(bcells, cells = cells_keep$barcode)

# HTO demux

bcells_demux <- HTODemux(bcells_pass, assay = "HTO", positive.quantile = 0.99)

bcells_pass_demux_meta <- bcells_demux@meta.data |>
    as_tibble(rownames = "barcode") |>
    filter(HTO_classification.global == "Singlet") |>
    inner_join(demuxlet_singlets, by = "barcode")

bcells_singlet <- subset(bcells_demux, cells = bcells_pass_demux_meta$barcode)

## PCA

bcells_singlet <- bcells_singlet |>
    FindVariableFeatures(nfeatures = 2000, selection.method = "vst") |>
    {function(x) ScaleData(x, features = rownames(x))}() |>
    {function(x) RunPCA(x, features = VariableFeatures(x))}() |>
    RunUMAP(dims = 1:20, verbose = FALSE) |>
    FindNeighbors(dims = 1:20, verbose = FALSE) |>
    FindClusters(resolution = 0.5, verbose = FALSE)

singlets_meta_data <- bcells_singlet@meta.data |>
    as_tibble(rownames = "barcode") |>
    select(barcode, cluster = seurat_clusters, stim = HTO_maxID, 
	   percent_mt, percent_ribo, n_genes = nFeature_RNA, umi = nCount_RNA) |>
    mutate(stim = factor(stim, levels = stims))

umap_df <- as.data.frame(bcells_singlet@reductions$umap@cell.embeddings) |>
    as_tibble(rownames = "barcode") |>
    left_join(singlets_meta_data)

umap_stims <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = stim)) +
    geom_point(size = .5) +
    scale_color_manual(values = stim_colors) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          legend.key.height = unit(.75, "lines")) +
    guides(color = guide_legend(override.aes = list(size = 2)))

ggsave("./plots/umap.png")
