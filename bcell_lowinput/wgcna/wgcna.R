# Data handling
library(tidyverse)
library(glue)

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


# Set colors
npg_colors <- ggsci::pal_npg()(10)

recolor <- function(x) { 
    case_when(x == "turquoise" ~ npg_colors[2],
	      x == "blue" ~ npg_colors[4],
	      x == "brown" ~ npg_colors[9],
	      x == "green" ~ npg_colors[3],
	      x == "red" ~ npg_colors[8],
	      x == "yellow" ~ "goldenrod1",
	      TRUE ~ x)
}


# Expression data
dat <- read_rds("./data/gene_expression.rds")

sample_table <- 
    tibble(sample_id = colnames(dat$counts)) |>
    extract(sample_id, c("group"), "[^_]+_(.+)", remove = FALSE) |>
    column_to_rownames("sample_id")

edger_res <- read_tsv("../results/edger/results.tsv")

stim_i <- "DN2"

## Use genes that are significant in the time course
sig_genes <- edger_res |>
    filter(stim == stim_i, FDR < 0.05) |>
    pull(gene_id)
    
## Create DESeq2 object for normalization of expression values.
## ~1 indicates no design, since we will not perform differential gene expression analysis here.
## Select samples with library size > 2MM reads.
## Filter for genes with normalized counts > 10 in at least 3 samples.
## Normalize counts using the variance-stabilizing transformation (VST)
dds <- 
    DESeqDataSetFromTximport(dat, sample_table, ~1) |>
    estimateSizeFactors() |>
    {function(x) x[, colSums(counts(x)) > 2e6]}() |>
    {function(x) x[rownames(x) %in% sig_genes, ]}() |>
    {function(x) x[, grepl(glue("_Unstim_0hrs|_{stim_i}_"), colnames(x))]}() |>
    {function(x) x[rowSums(counts(x, normalized = TRUE) >= 10 ) >= 3, ]}() |>
    vst()

## Transpose expression matrix to use with WGCNA
count_matrix <- t(assay(dds))

# Run WGCNA
allowWGCNAThreads()

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

ggsave(glue("./plots/wgcna_{stim_i}_sft.png"), 
       rsq_p | meank_p, height = 3, width = 7)

## Build WGCNA network 
## we need to reassing the 'cor()' function to avoid a bug in WGCNA
cut_height <- 0.3
beta_power <- ifelse(is.na(sft$powerEstimate), 30L, sft$powerEstimate)

cor <- WGCNA::cor

network <- 
    blockwiseModules(count_matrix,
		     corType = "bicor",
		     networkType = "signed",
		     TOMType = "signed",
		     power = beta_power, 
		     maxBlockSize = ncol(count_matrix),
		     minModuleSize = 30,
		     mergeCutHeight = cut_height,
		     detectCutHeight = 0.96,
		     numericLabels = FALSE,
		     randomSeed = 1,
		     verbose = 3)

cor <- stats::cor

## Gene classification into modules
module_genes <- enframe(network$colors, "gene_id", "module")

valid_modules <- module_genes |>
    count(module, sort = T) |>
    filter(module != "grey") |>
    pull(module)


## Module eigengenes
## A module eigengene is like a hypothetical gene the represents the module
module_eigen <- as_tibble(network$MEs, rownames = "sample_name")

# Save data for analysis of module preservation across stims
#write_tsv(module_genes, glue("./data/wgcna_{stim_i}_modules.tsv"))
#write_tsv(module_eigen, glue("./data/wgcna_{stim_i}_eigen.tsv"))
#write_rds(count_matrix, glue("./data/wgcna_{stim_i}_counts.rds"))

## plot dendogram to visualize modules
png(glue("./plots/wgcna_{stim_i}_dendro.png"), 
    width = 10, height = 6, unit = "in", res = 200)
plotDendroAndColors(network$dendrograms[[1]], 
		    recolor(network$colors),
		    "merged",
		    dendroLabels = FALSE,
		    addGuide = TRUE,
		    hang = 0.03, 
		    guideHang = 0.05)
dev.off()

# plot eigengene tree
eigengene_tree <- 
    as.dist(1 - bicor(select(network$MEs, -MEgrey))) |>
    hclust(method = "average")

png(glue("./plots/wgcna_{stim_i}_metree.png"), 
    width = 8, height = 8, unit = "in", res = 300)
plot(eigengene_tree, main = "Clustering of module eigengenes", xlab = "", sub = "")
abline(h = cut_height, col = "red")
dev.off()

# plot TOM
#sim_tom <- 
#    TOMsimilarityFromExpr(count_matrix, 
#			  corType = "bicor",
#			  networkType = "signed",
#			  power = beta_power,
#			  nThreads = 4)
#	
#diss_tom <- 1 - sim_tom
#gene_tree <- hclust(as.dist(diss_tom), method = "average")
#diag(diss_tom) <- NA
#
#png(glue("./plots/wgcna_{stim_i}_tom.png"), 
#    unit = "in", height = 6, width = 6, res = 200)
#TOMplot(diss_tom^8, gene_tree, Colors = network$colors, col = heat.colors(256))
#dev.off()


## Plot trends
scaled_mean_counts <- 
    count_matrix |>
    scale() |>
    as_tibble(rownames = "sample_name") |>
    separate(sample_name, c("subject_id", "stim", "timep"), sep = "_") |>
    mutate(timep = sub("hrs$", "", timep),
	   timep = factor(as.integer(timep)),
	   stim = factor(stim, levels = c("Unstim", stim_i))) |>
    pivot_longer(-c(subject_id, stim, timep), names_to = "gene_id") |>
    group_by(stim, timep, gene_id) |>
    summarise(mean_scaled_counts = mean(value)) |>
    ungroup()

kme_df <- 
    signedKME(count_matrix, orderMEs(network$MEs)) |>
    as_tibble(rownames = "gene_id") |>
    select(-kMEgrey) |>
    pivot_longer(-gene_id, names_to = "module")

top_kme <- kme_df |>
    filter(value > .7) |>
    group_by(gene_id) |>
    slice_max(value) |>
    group_by(module) |>
    top_n(500, abs(value)) |>
    ungroup() |>
    arrange(module, desc(abs(value)))

top_kme_counts <- 
    top_kme |>
    left_join(scaled_mean_counts, join_by(gene_id)) |>
    mutate(module = str_remove(module, "^kME"),
	   module = factor(module, levels = valid_modules))

module_colors <- setNames(recolor(valid_modules), valid_modules)

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
	   module = factor(module, levels = valid_modules))

trends_plot <- 
    ggplot(data = top_kme_counts,
	   aes(x = timep, y = mean_scaled_counts)) +
    geom_line(aes(group = gene_id, color = module), 
	      alpha = .05) +
    geom_line(data = mod_summ, 
	      aes(x = timep, y = value, group = module)) +
    geom_vline(xintercept = unique(top_kme_counts$timep), 
	       linetype = 2, linewidth = .1) +
    scale_color_manual(values = module_colors) +
    facet_wrap(~module, scale = "free", nrow = 2) +
    theme_bw() +
    theme(panel.grid = element_blank(),
	  legend.position = "none") +
    labs(x = "Hours", y = "Average scaled expression")

ggsave(glue("./plots/wgcna_{stim_i}_trends.png"), trends_plot, 
	    height = 3.5, width = 7)
   
# For presentation:
#pres_order <- 
#    c("Quick down and up", "Quick down and slow up", "Quick up – Quick down", 
#      "Quick up – peak at 24 hrs", "Late activation", "Constant down")
#
#top_kme_counts_pres <- 
#    top_kme_counts |>
#    mutate(facet_lab = fct_recode(module,
#			       "Constant down" = "turquoise",
#			       "Quick up – Quick down" = "blue",
#			       "Quick up – peak at 24 hrs" = "brown",
#			       "Late activation" = "yellow",
#			       "Quick down and up" = "green",
#			       "Quick down and slow up" = "red"),
#	   facet_lab = factor(facet_lab, levels = pres_order))
#
#mod_summ_pres <- 
#    mod_summ |>
#    mutate(facet_lab = fct_recode(module,
#			       "Constant down" = "turquoise",
#			       "Quick up – Quick down" = "blue",
#			       "Quick up – peak at 24 hrs" = "brown",
#			       "Late activation" = "yellow",
#			       "Quick down and up" = "green",
#			       "Quick down and slow up" = "red"),
#	   facet_lab = factor(facet_lab, levels = pres_order))
#   
#trends_plot_pres <- 
#    ggplot(data = top_kme_counts_pres,
#	   aes(x = timep, y = mean_scaled_counts)) +
#    geom_line(aes(group = gene_id, color = module), 
#	      alpha = .05) +
#    geom_line(data = mod_summ_pres, 
#	      aes(x = timep, y = value, group = module)) +
#    geom_vline(xintercept = unique(top_kme_counts$timep), 
#	       linetype = 2, linewidth = .1) +
#    scale_color_manual(values = module_colors) +
#    facet_wrap(~facet_lab, scale = "free", nrow = 2) +
#    theme_bw() +
#    theme(panel.grid = element_blank(),
#	  legend.position = "none") +
#    labs(x = "Hours", y = "Average scaled expression")
#
#ggsave(glue("./plots/wgcna_{stim_i}_trends_pres.png"), trends_plot_pres, 
#	    height = 3.5, width = 7)
#
#

# GO
run_enrichment <- function(gene_list) {

    enrichGO(gene = gene_list,
	OrgDb = org.Hs.eg.db,
	ont = "BP",
	keyType = "ENSEMBL",
	pAdjustMethod = "fdr",
	pvalueCutoff = 0.05,
	qvalueCutoff = 0.05,
	readable = TRUE)
}

## Run GO
## it takes a few minutes
go_res <- 
    top_kme |>
    mutate(gene_id = str_remove(gene_id, "\\.\\d+$"),
	   module = str_remove(module, "^kME")) |>
    filter(module != "grey") |>
    {function(x) split(x, x$module)}() |>
    map("gene_id") |>
    map(run_enrichment)

# Use module classification for GO instead of correlation between genes and MEs
#go_res2 <- 
#    module_genes |>
#    mutate(gene_id = str_remove(gene_id, "\\.\\d+$")) |>
#    filter(module != "grey") |>
#    {function(x) split(x, x$module)}() |>
#    map("gene_id") |>
#    map(run_enrichment)
#
#go_res2 |>
#    map_df(as_tibble, .id = "module") |>
#    select(module, Description, gene_id = geneID) |>
#    separate_rows(gene_id, sep = "/") |>
#    filter(gene_id == "STAT1")

go_res_top <- go_res |>
    map_df(as_tibble, .id = "module") |>
    mutate(Description = str_trunc(Description, width = 36)) |>
    group_by(module, Description) |>
    slice_max(-log10(pvalue)) |>
    group_by(module) |>
    top_n(10, -log10(pvalue)) |>
    ungroup() |>
    mutate(module = factor(module, levels = valid_modules))

go_plot <- 
    ggplot(data = go_res_top, 
       aes(x = "1", 
	   y = reorder_within(Description, by = -log10(pvalue), within = module))) +
    geom_point(aes(fill = -log10(pvalue), size = Count), shape = 21, stroke = .25) +
    scale_y_reordered() +
    scale_fill_gradient(low = "beige", high = "tomato4") +
    facet_wrap(~module, nrow = 2, scale = "free_y", drop = FALSE) +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
	  axis.text.y = element_text(size = 12),
	  axis.title = element_blank(),
	  panel.grid = element_blank(),
	  legend.position = "top",
	  legend.text = element_text(size = 14),
	  plot.background = element_rect(fill = "white", color = "white"),
	  plot.margin = margin(r = 1, l = 1, unit = "cm"),
	  strip.text = element_text(size = 14),
	  strip.clip = "off") +
    labs(fill = expression("-log"["10"]("p"))) +
    guides(fill = guide_colorbar(barheight = .7, barwidth = 10),
	   size = guide_legend(override.aes = list(fill = "black"))) +
    coord_cartesian(clip = "off")

ggsave(glue("./plots/wgcna_{stim_i}_go.png"), 
       go_plot, width = 11, height = 7, dpi = 300)
#

# Plot correlation between Lupus-associated genes and modules in a heatmap
# List of Lupus-associated genes
sle_genes <- 
    c("PTPN22", "FCGR2A", "TNFSF4", "NCF2", "SPRED2", "IFIH1", "STAT1", "STAT4",
      "IKZF1", "IKZF2", "IKZF3", "PXK", "IL12A", "BANK1", "BLK", "TCF7", "SKP1",
      "TNIP1", "MIR3142HG", "C4A", "C4B", "PRDM1", "TNFAIP3", "IRF5", "IRF7", "IRF8",
      "WDFY4", "ARID5B", "CD44", "ETS1", "SH2B3", "SLC15A4", "CSK", "CIITA", "SOCS1", 
      "CLEC16A", "ITGAM", "ITGAX", "GSDMB", "ORMDL3", "TYK2", "UBE2L3", "NR4A1", "MYC", "TREX1",
      "TLR7", "IRAK1")

gene_names <- distinct(edger_res, gene_id, gene_name)  

sle_genes_cormatrix <- 
    gene_names |>
    filter(gene_name %in% sle_genes) |>
    left_join(kme_df, join_by(gene_id)) |>
    mutate(module = str_remove(module, "^kME")) |>
    filter(!is.na(module), module != "grey") |>
    mutate(module = factor(module, levels = valid_modules)) |>
    select(-gene_id) |>
    pivot_wider(names_from = module, values_from = value) |>
    column_to_rownames("gene_name") |>
    data.matrix()

png(glue("./plots/wgcna_{stim_i}_slegenes.png"), 
    units = "in", height = 8, width = 4, res = 400)
pheatmap(sle_genes_cormatrix,
	 fontsize = 9,
	 color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
	 legend_breaks = seq(-1, 1, by = .25))
dev.off()

## Interferon genes
#pathways <- 
#    "/lab-share/IM-Gutierrez-e2/Public/References/msigdb/h.all.v2022.1.Hs.symbols.gmt.txt" |>
#    fgsea::gmtPathways() |>
#    {function(x) keep(x, grepl("interferon", names(x), ignore.case = TRUE))}()
#
#ifn_alpha_cormatrix <- 
#    gene_names |>
#    filter(gene_name %in% pathways[[1]]) |>
#    left_join(kme_df, join_by(gene_id)) |>
#    mutate(module = str_remove(module, "^kME")) |>
#    filter(!is.na(module), module != "grey") |>
#    select(-gene_id) |>
#    pivot_wider(names_from = module, values_from = value) |>
#    column_to_rownames("gene_name") |>
#    data.matrix()
#
#ifn_gamma_cormatrix <- 
#    gene_names |>
#    filter(gene_name %in% pathways[[2]]) |>
#    left_join(kme_df, join_by(gene_id)) |>
#    mutate(module = str_remove(module, "^kME")) |>
#    filter(!is.na(module), module != "grey") |>
#    select(-gene_id) |>
#    pivot_wider(names_from = module, values_from = value) |>
#    column_to_rownames("gene_name") |>
#    data.matrix()
#
#png(glue("./plots/wgcna_{stim_i}_ifna.png"), 
#    units = "in", height = 14, width = 4, res = 400)
#pheatmap(ifn_alpha_cormatrix,
#	 fontsize = 9,
#	 color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
#	 legend_breaks = seq(-1, 1, by = .25))
#dev.off()
#
#png(glue("./plots/wgcna_{stim_i}_ifng.png"), units = "in", height = 18, width = 4, res = 400)
#pheatmap(ifn_gamma_cormatrix,
#	 fontsize = 9,
#	 color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
#	 legend_breaks = seq(-1, 1, by = .25))
#dev.off()
#
#
#
#
#
## Hub genes
#hub_genes <- 
#    chooseTopHubInEachModule(count_matrix, 
#			     network$colors, 
#			     power = beta_power, 
#			     type = "signed") |>
#    enframe("module", "gene_id") |>
#    left_join(distinct(edger_res, gene_id, gene_name), join_by(gene_id))
#
#gene_adj <-
#    map_df(setNames(unique(module_genes$module), unique(module_genes$module)),
#	   ~adjacency(count_matrix[, network$colors == .x],
#		      power = beta_power, 
#		      type = "signed") |>
#	   rowSums() |>
#	   enframe("gene_id", "sum_a") |>
#	   arrange(desc(sum_a)),
#	   .id = "module") |>
#    left_join(distinct(edger_res, gene_id, gene_name), join_by(gene_id)) |>
#    select(module, gene_id, gene_name, sum_a)
#   
#gene_adj |>
#    group_split(module) |>
#    map(~head(., 20))
#





## Try dynamic tree cut instead of a fixed height cut to call modules
#dissTOM <- 1 - TOMsimilarityFromExpr(vst_dn2, power = 24, networkType = "signed")
#geneTree <- hclust(as.dist(dissTOM), method = "average" )
#
#dynamic_mods <-
#    cutreeDynamic(dendro = geneTree, distM = dissTOM,
#		  deepSplit = 2, pamRespectsDendro = FALSE,
#		  minClusterSize = 30)
#
#dynamic_colors <- labels2colors(dynamic_mods)
#
## merge similar modules
#me_list <- moduleEigengenes(vst_dn2, colors = dynamic_colors)
#mes <- me_list$eigengenes
#me_diss <- 1 - cor(mes)
#me_tree <- hclust(as.dist(me_diss), method = "average")
#
#png("./plots/wgcna_metree_dn2.png", width = 8, height = 6, unit = "in", res = 300)
#plot(me_tree, main = "Clustering of module eigengenes", xlab = "", sub = "")
#abline(h=0.25, col = "red")
#dev.off()
#
#
#mod_merge <- mergeCloseModules(vst_dn2, dynamic_colors, cutHeight = 0.25, verbose = 3)
#merged_colors <- mod_merge$colors
#merged_mes = mod_merge$newMEs
#
#png("./plots/wgcna_dendro_dynamic_dn2.png", width = 8, height = 4, unit = "in", res = 300)
#plotDendroAndColors(geneTree, cbind(dynamic_colors, merged_colors),
#                  c("Dynamic Tree Cut", "Merged dynamic"),
#                  dendroLabels = FALSE, hang = 0.03,
#                  addGuide = TRUE, guideHang = 0.05)
#dev.off()


## Network visualization
#green_genes <- module_genes |>
#    filter(module == "green") |>
#    pull(gene_id)
#
#adj_green <- adjacency(count_matrix[, green_genes], power = beta_power, type = "signed")
#tom_green <- TOMsimilarity(adj_green, TOMType = "signed")
#
#
#tom_green[1:5, 1:5]
#
#
#
