# ==============================================================================
# Description:  Performs differential expression analysis on the CITE-seq data 
#               using a pseudobulk approach. It aggregates single-cell counts 
#               (both RNA and ADT) at the donor/condition level, allowing us to 
#               use robust bulk RNA-seq statistical models (edgeR) while strictly 
#               accounting for paired biological donor variance.
# Input:        1. v4_seurat_qced.rds (Clean Seurat object)
#               2. features.tsv.gz (Gene annotations)
# Output:       Supplementary_Data_5_DGE_citeseq.xlsx
#               (Master Excel file containing ADT, RNA, and Cluster DGE results)
# ==============================================================================

# ==============================================================================
# 1. Environment Setup & Data Import
# ==============================================================================
library(Seurat)
library(tidyverse)
library(edgeR)
library(writexl)

# Load the finalized, QC-filtered Seurat object
bcells <- read_rds("../04_citeseq/2-processing/data/v4_seurat_qced.rds")

# Load gene annotations from the Cell Ranger output
features <- 
    "../04_citeseq/data/cellranger/1984/filtered_feature_bc_matrix/features.tsv.gz" |>
    read_tsv(col_names = c("gene_id", "gene_name", "phenotype"))

# Standardize the order of stimulation conditions for factor levels
fct_levels <- 
    c("Unstim", "IL4c24", "IL4c72", "TLR7c24", "TLR7c72",
      "BCRc24", "BCRc72", "DN2c24", "DN2c72")

# ==============================================================================
# 2. ADT (Surface Protein) Differential Expression
# ==============================================================================

# Aggregate single-cell ADT counts by human donor and stimulation condition.
# This yields a traditional bulk-like count matrix.
pseudo_adt <- 
    AggregateExpression(bcells,
			group.by = c("donor_id", "dmm_hto_call"),
			assays = "ADT",
			slot = "counts",
			return.seurat = FALSE)

pseudo_adt_mat <- pseudo_adt$ADT

# Parse the aggregated column names to build the edgeR design matrix
pseudo_adt_meta <-
    tibble(sample_id = colnames(pseudo_adt_mat)) |>
    separate(sample_id, c("donor_id", "condition"), sep = "_", remove = FALSE) |>
    mutate(donor_id = factor(donor_id),
	   condition = case_match(condition, 
			    "Unstim 0h" ~ "Unstim",
			    "IL4 24h" ~ "IL4c24",
			    "IL4 72h" ~ "IL4c72",
			    "TLR7 24h" ~ "TLR7c24",
			    "TLR7 72h" ~ "TLR7c72",
			    "BCR 24h" ~ "BCRc24",
			    "BCR 72h" ~ "BCRc72",
			    "DN2 24h" ~ "DN2c24",
			    "DN2 72h" ~ "DN2c72"),
	   condition = factor(condition, levels = fct_levels)
    )

# Initialize edgeR, calculate TMM normalization factors, and estimate dispersion.
# The design formula strictly pairs samples by donor_id to isolate condition effects.
y_adt <- DGEList(counts = pseudo_adt_mat, samples = pseudo_adt_meta) 
y_adt <- calcNormFactors(y_adt)
design_adt <- model.matrix(~ 0 + condition + donor_id, data = y_adt$samples)
colnames(design_adt) <- str_remove(colnames(design_adt), "condition")
y_adt <- estimateDisp(y_adt, design_adt)
fit_adt <- glmQLFit(y_adt, design_adt)

# Define pairwise contrasts of interest
my_contrasts_adt <- 
    makeContrasts(
		  IL4c24_vs_Unstim = IL4c24 - Unstim,
		  TLR7c24_vs_Unstim = TLR7c24 - Unstim,
		  BCRc24_vs_Unstim = BCRc24 - Unstim,
		  DN2c24_vs_Unstim = DN2c24 - Unstim,
		  
		  IL4c72_vs_Unstim = IL4c72 - Unstim,
		  TLR7c72_vs_Unstim = TLR7c72 - Unstim,
		  BCRc72_vs_Unstim = BCRc72 - Unstim,
		  DN2c72_vs_Unstim = DN2c72 - Unstim,

		  TLR7c24_vs_IL4c24 = TLR7c24 - IL4c24,
		  BCRc24_vs_IL4c24 = BCRc24 - IL4c24,
		  DN2c24_vs_IL4c24 = DN2c24 - IL4c24,
		  BCRc24_vs_TLR7c24 = BCRc24 - TLR7c24,
		  DN2c24_vs_TLR7c24 = DN2c24 - TLR7c24,

		  TLR7c72_vs_IL4c72 = TLR7c72 - IL4c72,
		  BCRc72_vs_IL4c72 = BCRc72 - IL4c72,
		  DN2c72_vs_IL4c72 = DN2c72 - IL4c72,
		  BCRc72_vs_TLR7c72 = BCRc72 - TLR7c72,
		  DN2c72_vs_TLR7c72 = DN2c72 - TLR7c72,

		  levels = design_adt
    )

# Run Quasi-Likelihood F-tests for all ADT contrasts and compile results
adt_res <- 
    colnames(my_contrasts_adt) |>
    set_names() |>
    map_dfr(function(i) glmQLFTest(fit_adt, contrast = my_contrasts_adt[, i]) |>
	    topTags(n = Inf) |>
	    {function(x) x$table}() |>
	    as_tibble(rownames = "gene_id"),
	    .id = "contrast") |>
    separate(contrast, c("stim1", "stim2"), sep = "_vs_") |>
    mutate_at(vars(stim1, stim2), ~factor(., levels = fct_levels)) |>
    arrange(stim1, stim2, gene_id) |>
    select(stim1, stim2, adt = gene_id, everything())

# ==============================================================================
# 3. RNA (Gene Expression) Differential Expression by Condition
# ==============================================================================

# Aggregate single-cell RNA counts by human donor and stimulation condition
pseudo_rna <- 
    AggregateExpression(bcells,
			group.by = c("donor_id", "dmm_hto_call"),
			assays = "RNA",
			slot = "counts",
			return.seurat = FALSE)

pseudo_rna_mat <- pseudo_rna$RNA

# Build metadata table identically to the ADT pipeline
pseudo_rna_meta <-
    tibble(sample_id = colnames(pseudo_rna_mat)) |>
    separate(sample_id, c("donor_id", "condition"), sep = "_", remove = FALSE) |>
    mutate(donor_id = factor(donor_id),
	   condition = case_match(condition, 
			    "Unstim 0h" ~ "Unstim",
			    "IL4 24h" ~ "IL4c24",
			    "IL4 72h" ~ "IL4c72",
			    "TLR7 24h" ~ "TLR7c24",
			    "TLR7 72h" ~ "TLR7c72",
			    "BCR 24h" ~ "BCRc24",
			    "BCR 72h" ~ "BCRc72",
			    "DN2 24h" ~ "DN2c24",
			    "DN2 72h" ~ "DN2c72"),
	   condition = factor(condition, levels = fct_levels) 
    )

# edgeR Modeling: Includes filterByExpr to remove lowly expressed genes across pseudobulks
y_rna <- DGEList(counts = pseudo_rna_mat, samples = pseudo_rna_meta) 
keep_rna <- filterByExpr(y_rna, group = y_rna$samples$condition)
y_rna <- y_rna[keep_rna, , keep.lib.sizes = FALSE]
y_rna <- calcNormFactors(y_rna)
design_rna <- model.matrix(~ 0 + condition + donor_id, data = y_rna$samples)
colnames(design_rna) <- str_remove(colnames(design_rna), "condition")
y_rna <- estimateDisp(y_rna, design_rna)
fit_rna <- glmQLFit(y_rna, design_rna)

# Define pairwise contrasts of interest
my_contrasts_rna <- 
    makeContrasts(
		  IL4c24_vs_Unstim = IL4c24 - Unstim,
		  TLR7c24_vs_Unstim = TLR7c24 - Unstim,
		  BCRc24_vs_Unstim = BCRc24 - Unstim,
		  DN2c24_vs_Unstim = DN2c24 - Unstim,
		  
		  IL4c72_vs_Unstim = IL4c72 - Unstim,
		  TLR7c72_vs_Unstim = TLR7c72 - Unstim,
		  BCRc72_vs_Unstim = BCRc72 - Unstim,
		  DN2c72_vs_Unstim = DN2c72 - Unstim,

		  TLR7c24_vs_IL4c24 = TLR7c24 - IL4c24,
		  BCRc24_vs_IL4c24 = BCRc24 - IL4c24,
		  DN2c24_vs_IL4c24 = DN2c24 - IL4c24,
		  BCRc24_vs_TLR7c24 = BCRc24 - TLR7c24,
		  DN2c24_vs_TLR7c24 = DN2c24 - TLR7c24,

		  TLR7c72_vs_IL4c72 = TLR7c72 - IL4c72,
		  BCRc72_vs_IL4c72 = BCRc72 - IL4c72,
		  DN2c72_vs_IL4c72 = DN2c72 - IL4c72,
		  BCRc72_vs_TLR7c72 = BCRc72 - TLR7c72,
		  DN2c72_vs_TLR7c72 = DN2c72 - TLR7c72,

		  levels = design_rna
    )

# Run QLF tests and merge with gene annotations
rna_res <- 
    colnames(my_contrasts_rna) |>
    set_names() |>
    map_dfr(function(i) glmQLFTest(fit_rna, contrast = my_contrasts_rna[, i]) |>
	    topTags(n = Inf) |>
	    {function(x) x$table}() |>
	    as_tibble(rownames = "gene_id"),
	    .id = "contrast") |>
    left_join(features, join_by(gene_id)) |>
    separate(contrast, c("stim1", "stim2"), sep = "_vs_") |>
    mutate_at(vars(stim1, stim2), ~factor(., levels = fct_levels)) |>
    arrange(stim1, stim2, gene_name) |>
    select(stim1, stim2, gene_id, gene_name, logFC:FDR)

# ==============================================================================
# 4. RNA Differential Expression by UMAP Cluster
# ==============================================================================

# Aggregate single-cell RNA counts by human donor and UMAP cluster identity
pseudo_cluster <- 
    AggregateExpression(bcells,
			group.by = c("donor_id", "seurat_clusters"),
			assays = "RNA",
			slot = "counts",
			return.seurat = FALSE)

pseudo_cluster_mat <- pseudo_cluster$RNA

pseudo_cluster_meta <-
    tibble(sample_id = colnames(pseudo_cluster_mat)) |>
    separate(sample_id, c("donor_id", "cluster"), sep = "_", remove = FALSE) |>
    mutate_at(vars(donor_id, cluster), factor)

# edgeR Modeling
y_cluster <- DGEList(counts = pseudo_cluster_mat, samples = pseudo_cluster_meta) 
keep_cluster <- filterByExpr(y_cluster, group = y_cluster$samples$cluster)
y_cluster <- y_cluster[keep_cluster, , keep.lib.sizes = FALSE]
y_cluster <- calcNormFactors(y_cluster)
design_cluster <- model.matrix(~ 0 + cluster + donor_id, data = y_cluster$samples)
y_cluster <- estimateDisp(y_cluster, design_cluster)
fit_cluster <- glmQLFit(y_cluster, design_cluster)

# Define cluster-specific contrasts.
# Tests expression in the target cluster against the unweighted mean expression 
# of all 13 other clusters combined.
cluster_contrasts <- 
    makeContrasts(
		  cluster0 = cluster0 - (cluster1 + cluster2 + cluster3 + cluster4 + cluster5 + cluster6 + cluster7 + cluster8 + cluster9 + cluster10 + cluster11 + cluster12 + cluster13) / 13,
		  cluster1 = cluster1 - (cluster0 + cluster2 + cluster3 + cluster4 + cluster5 + cluster6 + cluster7 + cluster8 + cluster9 + cluster10 + cluster11 + cluster12 + cluster13) / 13,
		  cluster2 = cluster2 - (cluster0 + cluster1 + cluster3 + cluster4 + cluster5 + cluster6 + cluster7 + cluster8 + cluster9 + cluster10 + cluster11 + cluster12 + cluster13) / 13,
		  cluster3 = cluster3 - (cluster0 + cluster1 + cluster2 + cluster4 + cluster5 + cluster6 + cluster7 + cluster8 + cluster9 + cluster10 + cluster11 + cluster12 + cluster13) / 13,
		  cluster4 = cluster4 - (cluster0 + cluster1 + cluster2 + cluster3 + cluster5 + cluster6 + cluster7 + cluster8 + cluster9 + cluster10 + cluster11 + cluster12 + cluster13) / 13,
		  cluster5 = cluster5 - (cluster0 + cluster1 + cluster2 + cluster3 + cluster4 + cluster6 + cluster7 + cluster8 + cluster9 + cluster10 + cluster11 + cluster12 + cluster13) / 13,
		  cluster6 = cluster6 - (cluster0 + cluster1 + cluster2 + cluster3 + cluster4 + cluster5 + cluster7 + cluster8 + cluster9 + cluster10 + cluster11 + cluster12 + cluster13) / 13,
		  cluster7 = cluster7 - (cluster0 + cluster1 + cluster2 + cluster3 + cluster4 + cluster5 + cluster6 + cluster8 + cluster9 + cluster10 + cluster11 + cluster12 + cluster13) / 13,
		  cluster8 = cluster8 - (cluster0 + cluster1 + cluster2 + cluster3 + cluster4 + cluster5 + cluster6 + cluster7 + cluster9 + cluster10 + cluster11 + cluster12 + cluster13) / 13,
		  cluster9 = cluster9 - (cluster0 + cluster1 + cluster2 + cluster3 + cluster4 + cluster5 + cluster6 + cluster7 + cluster8 + cluster10 + cluster11 + cluster12 + cluster13) / 13,
		  cluster10 = cluster10 - (cluster0 + cluster1 + cluster2 + cluster3 + cluster4 + cluster5 + cluster6 + cluster7 + cluster8 + cluster9 + cluster11 + cluster12 + cluster13) / 13,
		  cluster11 = cluster11 - (cluster0 + cluster1 + cluster2 + cluster3 + cluster4 + cluster5 + cluster6 + cluster7 + cluster8 + cluster9 + cluster10 + cluster12 + cluster13) / 13,
		  cluster12 = cluster12 - (cluster0 + cluster1 + cluster2 + cluster3 + cluster4 + cluster5 + cluster6 + cluster7 + cluster8 + cluster9 + cluster10 + cluster11 + cluster13) / 13,
		  cluster13 = cluster13 - (cluster0 + cluster1 + cluster2 + cluster3 + cluster4 + cluster5 + cluster6 + cluster7 + cluster8 + cluster9 + cluster10 + cluster11 + cluster12) / 13,
		  levels = design_cluster
    )

# Run QLF tests and compile results
cluster_res <-     
    colnames(cluster_contrasts) |>
    set_names() |>
    map_dfr(function(i) glmQLFTest(fit_cluster, contrast = cluster_contrasts[, i]) |>
	    topTags(n = Inf) |>
	    {function(x) x$table}() |>
	    as_tibble(rownames = "gene_id"),
	    .id = "cluster") |>
    left_join(features, join_by(gene_id)) |>
    select(cluster, gene_id, gene_name, logFC:FDR)

# ==============================================================================
# 5. Data Export
# ==============================================================================
supp_data_igd_igm <- 
    list("ADT_condition" = adt_res, 
	 "RNA_condition" = rna_res, 
	 "RNA_cluster" = cluster_res)

write_xlsx(supp_data_igd_igm, "./data/Supplementary_Data_5_DGE_citeseq.xlsx")
