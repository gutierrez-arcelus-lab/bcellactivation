# Packages
## single-cell data analysis
library(Seurat)

## Data wrangling
library(tidyverse)

## Plotting
library(tidytext)
library(RColorBrewer)
library(scico)
library(ggsci)
library(cowplot)
library(sclibr)


# Directories
labshr <- "/lab-share/IM-Gutierrez-e2/Public/Lab_datasets"

pilot1_dir <- 
    file.path(labshr, "CITEseq_pilot", 
              "SN0231064/KW9100_Maria/210726_10X_KW9100-2_bcl/cellranger-6.0.1",
              "GRCh38/BRI-1283/outs/filtered_feature_bc_matrix")

pilot2_dir <- 
    file.path(labshr, "CITEseq_pilot_2",
	      "SN0257788/broad/hptmp/curtism/bwh10x/KW10170_Maria/220617_10X_KW10170_bcl",
	      "cellranger-6.1.1/GRCh38/BRI-1743_hashing/outs/filtered_feature_bc_matrix")

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
pilot1_features <- file.path(pilot1_dir, "features.tsv.gz") |>
    read_tsv(col_names = c("gene_id", "gene_name", "phenotype"))

pilot2_features <- file.path(pilot2_dir, "features.tsv.gz") |>
    read_tsv(col_names = c("gene_id", "gene_name", "phenotype"))

lib1984_features <- file.path(lib1984_dir, "features.tsv.gz") |>
    read_tsv(col_names = c("gene_id", "gene_name", "phenotype"))

lib1988_features <- file.path(lib1988_dir, "features.tsv.gz") |>
    read_tsv(col_names = c("gene_id", "gene_name", "phenotype"))

lib1990_features <- file.path(lib1990_dir, "features.tsv.gz") |>
    read_tsv(col_names = c("gene_id", "gene_name", "phenotype"))


# Mitochondrial and Ribosomal protein genes
# Gene IDs are the same across CITE-seq runs, so we'll use pilot 1.
mt_genes <- pilot1_features |>
    filter(grepl("^MT-", gene_name)) |>
    pull(gene_id)

ribo_genes <- pilot1_features |>
    filter(grepl("^RPS\\d+|^RPL\\d+", gene_name)) |>
    pull(gene_id)

# Hashtag IDs
pilot1_stims <- 
    c("Hashtag1" = "BCR 72hr",
      "Hashtag2" = "TLR7 72hr",
      "Hashtag3" = "BCR 24hr", 
      "Hashtag4" = "TLR7 24hr",
      "Hashtag5" = "Res 24hr",
      "Hashtag6" = "Day 0")

pilot2_stims <- 
    c("Hashtag1" = "Day 0", 
      "Hashtag2" = "IL4 24hr",
      "Hashtag3" = "BCR 24hr",
      "Hashtag4" = "BCR+TLR7 24hr",
      "Hashtag5" = "TLR7 24hr", 
      "Hashtag6" = "sCD40L 24hr",
      "Hashtag7" = "CpG 24hr",
      "Hashtag8" = "DN2 24hr",
      "Hashtag9" = "BCR 72hr",
      "Hashtag10" = "BCR+TLR7 72hr",
      "Hashtag12" = "TLR7 72hr",
      "Hashtag13" = "sCD40L 72hr",
      "Hashtag14" = "DN2 72hr",
      "Hashtag15" = "CpG 72hr")

mgb_stims <- 
    c("Hashtag6" = "Day 0",
      "Hashtag7" = "IL4 24h",
      "Hashtag8" = "IL4 72h",
      "Hashtag9" = "BCR 24h",
      "Hashtag10" = "BCR 72h",
      "Hashtag12" = "TLR7 24h",
      "Hashtag13" = "TLR7 72h",
      "Hashtag14" = "DN2 24h",
      "Hashtag15" = "DN2 72h")


# Seurat objects
# using function from custom R package 'sclibr'

pilot1_obj <- make_seurat(pilot1_dir, project_id = "pilot1", hto_names = pilot1_stims, 
			  mito_ids = mt_genes, ribo_ids = ribo_genes)

pilot2_obj <- make_seurat(pilot2_dir, project_id = "pilot2", hto_names = pilot2_stims, 
			  mito_ids = mt_genes, ribo_ids = ribo_genes)

lib1984_obj <- make_seurat(lib1984_dir, project_id = "1984", hto_names = mgb_stims,
			   mito_ids = mt_genes, ribo_ids = ribo_genes)

lib1988_obj <- make_seurat(lib1988_dir, project_id = "1988", hto_names = mgb_stims,
			   mito_ids = mt_genes, ribo_ids = ribo_genes)

lib1990_obj <- make_seurat(lib1990_dir, project_id = "1990", hto_names = mgb_stims,
			   mito_ids = mt_genes, ribo_ids = ribo_genes)


