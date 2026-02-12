# Data handling
library(tidyverse)
library(glue)
library(furrr)

# import data
library(tximport)
library(DESeq2)

# Run main analyses
library(WGCNA)

# Gene Ontology analysis
library(clusterProfiler)
library(org.Hs.eg.db)

# Plotting
library(pheatmap)
library(RColorBrewer)
library(tidytext)
library(patchwork)

count <- dplyr::count
select <- dplyr::select

# Set WGCNA cores
allowWGCNAThreads(parallelly::availableCores())

# Expression data
dat <- read_rds("./data/gene_expression.rds")

sample_table <- 
    tibble(sample_id = colnames(dat$counts)) |>
    extract(sample_id, c("group"), "[^_]+_(.+)", remove = FALSE) |>
    column_to_rownames("sample_id")

edger_res <- read_tsv("../results/edger/results.tsv")
gene_names <- distinct(edger_res, gene_id, gene_name)  

## Create DESeq2 object for normalization of expression values.
## ~1 indicates no design, since we will not perform differential gene expression analysis here.
## Normalize counts using the variance-stabilizing transformation (VST)
dds <- 
    DESeqDataSetFromTximport(dat, sample_table, ~1) |>
    estimateSizeFactors() |>
    {function(x) x[rownames(x) %in% unique(edger_res$gene_id), ]}() |>
    {function(x) x[, grepl("_Unstim_0hrs|_CD40L_|_TLR7_|_TLR9_|_BCR_|_BCR-TLR7_|_DN2_", colnames(x))]}() |>
    vst()

## Transpose expression matrix to use with WGCNA
count_matrix <- t(assay(dds))

# Run WGCNA

## Explore multiple values of Beta (power) to find the one that satisfies scale-free topology
betas <- c(1:10, seq(12, 30, by = 2))

sft <- 
    pickSoftThreshold(count_matrix, 
		      powerVector = betas, 
		      networkType = "signed", 
		      blockSize = ncol(count_matrix),
		      verbose = 3)

## Plot beta vs r^2 between connectivity and frequency of connectivity
## and beta vs mean connectivity.
## The goal here is to find the minimum value of beta that maximizes r^2
## and minimizes connectivity.
## That's because we believe that, in gene networks, there should be few genes
## with many connections (hub genes), and most genes should connect to a small
## number of genes.
rsq_p <- 
    ggplot(data = sft$fitIndices, 
	   aes(x = Power, y = SFT.R.sq)) +
    geom_hline(yintercept = 0.85, color = "red") +
    geom_point() +
    geom_line() +
    scale_y_continuous(limits = c(0, 1),
		       breaks = seq(0, 1, by = .2)) +
    theme_bw()

meank_p <-
    ggplot(data = sft$fitIndices,
	   aes(x = Power, y = mean.k.)) +
    geom_point() +
    geom_line() +
    theme_bw()

ggsave(glue("./plots/sft_all.png"), 
       rsq_p | meank_p, height = 3, width = 7)

### Build WGCNA network 
### we need to reassing the 'cor()' function to avoid a bug in WGCNA
#merge_cut_height <- 0.2
#detect_cut_height <- 0.99

beta_power <- ifelse(is.na(sft$powerEstimate), 30L, sft$powerEstimate)

cor <- WGCNA::cor

network <- 
    blockwiseModules(count_matrix,
		     corType = "pearson",
		     networkType = "signed",
		     TOMType = "signed",
		     power = beta_power, 
		     maxBlockSize = ncol(count_matrix),
		     numericLabels = FALSE,
		     detectCutHeight = 0.99,
		     mergeCutHeight = .2,
		     minKMEtoStay = .5,
		     minModuleSize = 30,
		     randomSeed = 1,
		     verbose = 3)

cor <- stats::cor

module_genes <- enframe(network$colors, "gene_id", "module")

module_colors <- 
    module_genes |>
    count(module, sort = TRUE) |>
    filter(module != "grey") |>
    select(-n) |>
    deframe()

## plot dendogram to visualize modules
png(glue("./plots/dendro_all.png"), 
    width = 10, height = 6, unit = "in", res = 200)
plotDendroAndColors(network$dendrograms[[1]], 
		    network$colors,
		    "merged",
		    dendroLabels = FALSE,
		    addGuide = TRUE,
		    hang = 0.03, 
		    guideHang = 0.05,
		    abHeight = 0.99)
dev.off()

# plot eigengene tree
eigengene_tree <- 
    as.dist(1 - cor(select(network$MEs, -MEgrey))) |>
    hclust(method = "average")

png(glue("./plots/metree_all.png"), 
    width = 8, height = 8, unit = "in", res = 300)
plot(eigengene_tree, main = "Clustering of module eigengenes", xlab = "", sub = "")
abline(h = 0.2, col = "red")
dev.off()

# Gene classification into modules

## kME
kme_all_df <- 
    signedKME(count_matrix, 
	      orderMEs(network$MEs), 
	      outputColumnName = "") |>
    as_tibble(rownames = "gene_id") |>
    select(-grey) |>
    pivot_longer(-gene_id, names_to = "module", values_to = "kme")


## Module eigengenes
## A module eigengene is like a hypothetical gene the represents the module
module_eigen <- as_tibble(network$MEs, rownames = "sample_name")

## Plot trends
### Module assignment and higher kME are not always concordant.
### While module eigengene is based on original assigments (based on adjacency not correlation).
### Some genes have very low kME with the assigned module, and the threshold used in blockwideModules is not working properly.
### I'm selecting top 500 genes by highest kME with a module to test with GO
### By suggestion of WGCNA main author.
scaled_mean_counts <- 
    count_matrix |>
    scale() |>
    as_tibble(rownames = "sample_name") |>
    separate(sample_name, c("subject_id", "stim", "timep"), sep = "_") |>
    mutate(timep = sub("hrs$", "", timep),
	   timep = factor(as.integer(timep)),
	   stim = factor(stim, levels = c("Unstim", "CD40L", "TLR7", "TLR9", "BCR-TLR7", "BCR", "DN2"))) |>
    pivot_longer(-c(subject_id, stim, timep), names_to = "gene_id") |>
    group_by(timep, gene_id) |>
    summarise(mean_scaled_counts = mean(value)) |>
    ungroup()

kme_df <- 
    kme_all_df |>
    filter(kme >= 0.85) |>
    group_by(gene_id) |>
    slice_max(kme) |>
    ungroup()

top_kme <- kme_df |>
    group_by(module) |>
    top_n(500, kme) |>
    ungroup() |>
    arrange(module, desc(kme))

top_kme_counts <- 
    top_kme |>
    left_join(scaled_mean_counts, join_by(gene_id)) |>
    mutate(module = factor(module, levels = module_colors))

mod_summ <- module_eigen |>
    column_to_rownames("sample_name") |>
    scale() |>
    as_tibble(rownames = "sample_name") |>
    pivot_longer(-sample_name, names_to = "module") |>
    mutate(module = str_remove(module, "^ME")) |>
    separate(sample_name, c("donor_id", "stim", "timep"), sep = "_") |>
    group_by(timep, module) |>
    summarise(value = mean(value)) |>
    ungroup() |>
    filter(module != "grey") |>
    mutate(timep = str_remove(timep, "hrs$"),
	   timep = factor(as.integer(timep)),
	   module = factor(module, levels = module_colors))

trends_plot <- 
    ggplot(data = top_kme_counts,
	   aes(x = timep, y = mean_scaled_counts)) +
    geom_line(aes(group = gene_id, color = module), 
	      alpha = .05) +
    geom_line(data = mod_summ, 
	      aes(x = timep, y = value, group = module)) +
    geom_vline(xintercept = unique(top_kme_counts$timep), 
	       linetype = 2, linewidth = .1) +
    scale_color_manual(values = setNames(module_colors, module_colors)) +
    facet_wrap(~module, scale = "free", nrow = 2) +
    theme_bw() +
    theme(panel.grid = element_blank(),
	  legend.position = "none") +
    labs(x = "Hours", y = "Average scaled expression")

ggsave(glue("./plots/trends_all.png"), trends_plot, 
	    height = 3.5, width = 7)




traits <- 
    module_eigen |>
    select(sample_name) |> 
    separate(sample_name, c("donor_id", "stim", "timep"), sep = "_", remove = FALSE) |>
    mutate(i = 1) |>
    select(sample_name, stim, i) |>
    pivot_wider(names_from = stim, values_from = i) |> 
    select(sample_name, CD40L, TLR7, TLR9, BCR, `BCR-TLR7`, DN2) |>
    mutate_at(vars(-sample_name), ~replace_na(., 0))
    

module_trait_cor <-
    cor(select(module_eigen, -sample_name),
	select(traits, -sample_name), 
	use = "p")

module_trait_pvalue <- corPvalueStudent(module_trait_cor, nrow(module_eigen))

png(glue("./plots/module_trait_heatmap.png"), 
    width = 7, height = 7, unit = "in", res = 300)
par(mar = c(5, 8, 4, 2))
labeledHeatmap(Matrix = module_trait_cor,
               xLabels = colnames(module_trait_cor),
               yLabels = rownames(module_trait_cor),
               ySymbols = rownames(module_trait_cor) ,
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = paste(signif(module_trait_cor, 2), "\n(",
                                  signif(module_trait_pvalue, 1), ")", sep = ""),
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-Trait Relationships"))
dev.off()


## Intramodule connectivity
kim <-
    intramodularConnectivity.fromExpr(count_matrix, 
				      network$colors,
				      networkType = "signed",
				      power = beta_power,
				      scaleByMax = TRUE) 

kim_df <- 
    tibble(gene_id = colnames(count_matrix), kim = kim$kWithin) |>
    left_join(module_genes, join_by(gene_id)) |>
    left_join(gene_names) |>
    select(gene_id, gene_name, module, kim) |>
    filter(!is.na(kim))

