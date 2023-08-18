library(tidyverse)
library(tximport)
library(DESeq2)
library(WGCNA)
library(patchwork)

count <- dplyr::count
select <- dplyr::select

# Sample meta data
meta_data <- read_tsv("./data/sample_decode.tsv") |>
    separate(sample_name, c("name", "stim", "time"), sep = "_", remove = FALSE) |>
    mutate(time = factor(time, levels = str_sort(unique(time), numeric = TRUE))) |>
    arrange(stim, time, name)

sample_table <- meta_data |>
    unite("group", c(stim, time), sep = ".") |>
    select(sample_id, group) |>
    column_to_rownames("sample_id")

# Transcript to Gene map
tx_to_gene <- 
    file.path("/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens",
	      "gencode.v38.primary_assembly.annotation.gtf") |>
    read_tsv(comment = "#", col_names = FALSE) |>
    filter(X3 == "transcript") |>
    mutate(tx_id = str_extract(X9, "(?<=transcript_id\\s\")[^\"]+"),
	   gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
	   gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+")) |>
    select(tx_id, gene_id, gene_name)

# Expression data
salmon_files <- 
    sprintf("./results/salmon/%s/quant.sf", meta_data$sample_id) |>
    setNames(meta_data$sample_id)

txi <- tximport(salmon_files, 
		type = "salmon", 
		tx2gene = select(tx_to_gene, tx_id, gene_id))

dds <- DESeqDataSetFromTximport(txi, sample_table, ~1)
dds <- estimateSizeFactors(dds)
keep_y <- rowSums(counts(dds, normalized = TRUE) >= 10 ) >= 3
dds <- dds[keep_y, ]
dds <- collapseReplicates(dds, meta_data$sample_name, meta_data$sample_id)

gsg <- goodSamplesGenes(t(counts(dds)))
gsg$allOK

dds <- vst(dds)

# WGCNA
allowWGCNAThreads()

vst_counts <- t(assay(dds))

powers <- c(1:10, seq(12, 50, by = 2))

sft <- 
    pickSoftThreshold(vst_counts, 
		      powerVector = powers, 
		      networkType = "signed", 
		      blockSize = 30000,
		      verbose = 5)

rsq_p <- 
    ggplot(data = sft$fitIndices, 
	   aes(x = Power, y = SFT.R.sq)) +
    geom_point() +
    geom_line() +
    theme_bw()

meank_p <-
    ggplot(data = sft$fitIndices,
	   aes(x = Power, y = mean.k.)) +
    geom_point() +
    geom_line() +
    theme_bw()

ggsave("./plots/wgcna_powers.png", rsq_p | meank_p, height = 3)

cor <- WGCNA::cor

bwnet <- 
    blockwiseModules(vst_counts,
		     corType = "bicor",
		     networkType = "signed",
		     power = 14, 
		     maxBlockSize = 30000,
		     minModuleSize = 30,
		     mergeCutHeight = .25,
		     numericLabels = FALSE,
		     randomSeed = 1,
		     verbose = 3)

cor <- stats::cor

mod_eigengenes <- as_tibble(bwnet$MEs, rownames = "sample_name")

mod_gene <- enframe(bwnet$colors, "gene_id", "module")

mod_gene |> dplyr::count(module, sort = TRUE)

png("./plots/wgcna_dendrogram.png", width = 8, height = 4, unit = "in", res = 300)
plotDendroAndColors(bwnet$dendrograms[[1]], 
		    bwnet$colors,
		    "merged",
		    dendroLabels = FALSE,
		    addGuide = TRUE,
		    hang = 0.03, 
		    guideHang = 0.05)
dev.off()

mod_eigen_df <- mod_eigengenes |>
    pivot_longer(-sample_name, names_to = "module") |>
    separate(sample_name, c("subject_id", "stim", "timep"), sep = "_") |>
    mutate(timep = sub("hrs$", "", timep),
	   timep = factor(as.integer(timep)),
	   stim = factor(stim, levels = c("Unstim", "IL4", "CD40L", "BCR", "TLR7", "TLR9", "BCR-TLR7", "DN2"))) |>
    group_by(stim, timep, module) |>
    summarise(value = mean(value)) |>
    ungroup()

tp <- 
    ggplot(data = mod_eigen_df,
	   aes(x = timep, y = module)) +
    geom_tile(aes(fill = value)) +
    scale_fill_gradient2(high = "tomato3", mid = "white", low = "midnightblue",
			 guide = guide_colorbar(barheight = .5, barwidth = 10)) +
    facet_grid(~stim, scale = "free_x", space = "free") +
    theme_bw() +
    theme(legend.position = "top") +
    labs(x = "Hours", fill = NULL)

ggsave("./plots/wgcna_tp.png", tp, height = 3, width = 10)

count_norm_df <- vst_counts |>
    as_tibble(rownames = "sample_name") |>
    separate(sample_name, c("subject_id", "stim", "timep"), sep = "_") |>
    mutate(timep = sub("hrs$", "", timep),
	   timep = factor(as.integer(timep)),
	   stim = factor(stim, levels = c("Unstim", "IL4", "CD40L", "BCR", "TLR7", "TLR9", "BCR-TLR7", "DN2"))) |>
    pivot_longer(-c(subject_id, stim, timep), names_to = "gene_id") |>
    left_join(mod_gene, join_by(gene_id)) |>
    left_join(distinct(tx_to_gene, gene_id, gene_name), join_by(gene_id)) |>
    select(subject_id, stim, timep, gene_id, gene_name, module, value)

count_mod_summ <- count_norm_df |>
    group_by(stim, timep, gene_id, gene_name, module) |>
    summarise(value = mean(value)) |>
    ungroup()

yellow_p <- 
    ggplot(data = count_mod_summ |> filter(module == "yellow"),
       aes(x = timep, y = value)) +
    geom_line(aes(group = gene_id), alpha = .25) +
    facet_grid(~stim, scale = "free_x", space = "free") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(x = "Hours")

ggsave("./plots/wgcna_yellow.png", yellow_p, height = 3, width = 12)

black_p <- 
    ggplot(data = count_mod_summ |> filter(module == "black"),
       aes(x = timep, y = value)) +
    geom_line(aes(group = gene_id), alpha = .1) +
    facet_grid(~stim, scale = "free_x", space = "free") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(x = "Hours")

ggsave("./plots/wgcna_black.png", yellow_p, height = 3, width = 12)


count_mod_summ |> filter(gene_name == "TNFSF4") |> print(n = Inf)

# VST transform seems to bring very different expression levels close
# so samples don't look different anymore.
# Try TPM?
count_mod_summ |> filter(gene_name == "STAT1") |> print(n = Inf)


# Try TPM
keep_samples <- txi$counts |> 
    colSums() |>
    enframe("sample_id", "lib_size") |>
    left_join(meta_data, join_by(sample_id)) |>
    group_by(sample_name) |>
    slice_max(lib_size) |>
    ungroup()

tpm_matrix <- txi$counts[, keep_samples$sample_id]
colnames(tpm_matrix) <- keep_samples$sample_name
keep_tpm <- rowSums(tpm_matrix >= 10 ) >= 3
tpm_matrix <- tpm_matrix[keep_tpm, ]
tpm_matrix <- t(tpm_matrix)

powers <- c(1:10, seq(12, 50, by = 2))

sft <- 
    pickSoftThreshold(tpm_matrix, 
		      powerVector = powers, 
		      networkType = "signed", 
		      blockSize = 30000,
		      verbose = 5)

rsq_p <- 
    ggplot(data = sft$fitIndices, 
	   aes(x = Power, y = SFT.R.sq)) +
    geom_point() +
    geom_line() +
    theme_bw()

meank_p <-
    ggplot(data = sft$fitIndices,
	   aes(x = Power, y = mean.k.)) +
    geom_point() +
    geom_line() +
    theme_bw()

ggsave("./plots/wgcna_powers_tpm.png", rsq_p | meank_p, height = 3)

cor <- WGCNA::cor

bwnet_tpm <- 
    blockwiseModules(tpm_matrix,
		     corType = "bicor",
		     networkType = "signed",
		     power = 16, 
		     maxBlockSize = 30000,
		     minModuleSize = 30,
		     mergeCutHeight = .25,
		     numericLabels = FALSE,
		     randomSeed = 1,
		     verbose = 3)

cor <- stats::cor

mod_eigengenes_tpm <- as_tibble(bwnet_tpm$MEs, rownames = "sample_name")

mod_gene_tpm <- enframe(bwnet_tpm$colors, "gene_id", "module")

mod_gene_tpm |> dplyr::count(module, sort = TRUE)

png("./plots/wgcna_dendrogram_tpm.png", width = 8, height = 4, unit = "in", res = 300)
plotDendroAndColors(bwnet_tpm$dendrograms[[1]], 
		    bwnet_tpm$colors,
		    "merged",
		    dendroLabels = FALSE,
		    addGuide = TRUE,
		    hang = 0.03, 
		    guideHang = 0.05)
dev.off()

mod_eigen_tpm_df <- mod_eigengenes_tpm |>
    pivot_longer(-sample_name, names_to = "module") |>
    separate(sample_name, c("subject_id", "stim", "timep"), sep = "_") |>
    mutate(timep = sub("hrs$", "", timep),
	   timep = factor(as.integer(timep)),
	   stim = factor(stim, levels = c("Unstim", "IL4", "CD40L", "BCR", "TLR7", "TLR9", "BCR-TLR7", "DN2"))) |>
    group_by(stim, timep, module) |>
    summarise(value = mean(value)) |>
    ungroup()

tp_tpm <- 
    ggplot(data = mod_eigen_tpm_df,
	   aes(x = timep, y = module)) +
    geom_tile(aes(fill = value)) +
    scale_fill_gradient2(high = "tomato3", mid = "white", low = "midnightblue",
			 guide = guide_colorbar(barheight = .5, barwidth = 10)) +
    facet_grid(~stim, scale = "free_x", space = "free") +
    theme_bw() +
    theme(legend.position = "top") +
    labs(x = "Hours", fill = NULL)

ggsave("./plots/wgcna_tp_tpm.png", tp, height = 3, width = 10)


#DN2-only network 

## use genes that are significant in the time course
edger_res <- read_tsv("./results/edger/results.tsv")

dn2_sig <- edger_res |>
    filter(stim == "DN2", FDR < 0.05)

dds_dn2 <- DESeqDataSetFromTximport(txi, sample_table, ~1)
dds_dn2 <- estimateSizeFactors(dds_dn2)
dds_dn2 <- dds_dn2[dn2_sig$gene_id, ]
dds_dn2 <- collapseReplicates(dds_dn2, meta_data$sample_name, meta_data$sample_id)

dds_dn2 <- dds_dn2[, grepl("_DN2_|_Unstim_0hrs", colnames(dds_dn2))]

dds_dn2 <- vst(dds_dn2)
vst_dn2 <- t(assay(dds_dn2))

powers <- c(1:10, seq(12, 50, by = 2))

sft_dn2 <- 
    pickSoftThreshold(vst_dn2, 
		      powerVector = powers, 
		      networkType = "signed", 
		      blockSize = 30000,
		      verbose = 5)

rsq_dn2_p <- 
    ggplot(data = sft_dn2$fitIndices, 
	   aes(x = Power, y = SFT.R.sq)) +
    geom_point() +
    geom_line() +
    theme_bw()

meank_dn2_p <-
    ggplot(data = sft_dn2$fitIndices,
	   aes(x = Power, y = mean.k.)) +
    geom_point() +
    geom_line() +
    theme_bw()

ggsave("./plots/wgcna_powers_dn2.png", rsq_dn2_p | meank_dn2_p, height = 3)

cor <- WGCNA::cor

bwnet_dn2 <- 
    blockwiseModules(vst_dn2,
		     corType = "bicor",
		     networkType = "signed",
		     power = 24, 
		     maxBlockSize = 30000,
		     minModuleSize = 30,
		     mergeCutHeight = .25,
		     numericLabels = FALSE,
		     randomSeed = 1,
		     verbose = 3)

cor <- stats::cor

mod_eigengenes_dn2 <- as_tibble(bwnet_dn2$MEs, rownames = "sample_name")

mod_gene_dn2 <- enframe(bwnet_dn2$colors, "gene_id", "module")

mod_gene_dn2 |> dplyr::count(module, sort = TRUE)

png("./plots/wgcna_dendro_dn2.png", width = 8, height = 4, unit = "in", res = 300)
plotDendroAndColors(bwnet_dn2$dendrograms[[1]], 
		    bwnet_dn2$colors,
		    "merged",
		    dendroLabels = FALSE,
		    addGuide = TRUE,
		    hang = 0.03, 
		    guideHang = 0.05)
dev.off()

mod_eigen_dn2_df <- mod_eigengenes_dn2 |>
    pivot_longer(-sample_name, names_to = "module") |>
    separate(sample_name, c("subject_id", "stim", "timep"), sep = "_") |>
    mutate(timep = sub("hrs$", "", timep),
	   timep = factor(as.integer(timep)),
	   stim = factor(stim, levels = c("Unstim", "IL4", "CD40L", "BCR", "TLR7", "TLR9", "BCR-TLR7", "DN2"))) |>
    group_by(stim, timep, module) |>
    summarise(value = mean(value)) |>
    ungroup()

tp_dn2 <- 
    ggplot(data = mod_eigen_dn2_df,
	   aes(x = timep, y = module)) +
    geom_tile(aes(fill = value)) +
    scale_fill_gradient2(high = "tomato3", mid = "white", low = "midnightblue",
			 guide = guide_colorbar(barheight = .5, barwidth = 10)) +
    facet_grid(~stim, scale = "free_x", space = "free") +
    theme_bw() +
    theme(legend.position = "top") +
    labs(x = "Hours", fill = NULL)

ggsave("./plots/wgcna_tp_dn2.png", tp_dn2, height = 3, width = 10)

count_norm_dn2_df <- vst_dn2 |>
    as_tibble(rownames = "sample_name") |>
    separate(sample_name, c("subject_id", "stim", "timep"), sep = "_") |>
    mutate(timep = sub("hrs$", "", timep),
	   timep = factor(as.integer(timep)),
	   stim = factor(stim, levels = c("Unstim", "IL4", "CD40L", "BCR", "TLR7", "TLR9", "BCR-TLR7", "DN2"))) |>
    pivot_longer(-c(subject_id, stim, timep), names_to = "gene_id") |>
    left_join(mod_gene_dn2, join_by(gene_id)) |>
    left_join(distinct(tx_to_gene, gene_id, gene_name), join_by(gene_id)) |>
    select(subject_id, stim, timep, gene_id, gene_name, module, value)

count_mod_summ_dn2 <- count_norm_dn2_df |>
    group_by(stim, timep, gene_id, gene_name, module) |>
    summarise(value = mean(value)) |>
    ungroup()

mod_dn2_p <- 
    ggplot(data = count_mod_summ_dn2,
       aes(x = timep, y = value)) +
    geom_line(aes(group = gene_id), alpha = .1) +
    facet_grid(module~.) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(x = "Hours")

ggsave("./plots/wgcna_dn2.png", mod_dn2_p, height = 10, width = 5)


count_mod_summ_dn2 |> filter(gene_name == "STAT1")
count_mod_summ_dn2 |> filter(gene_name == "PTPN22")
count_mod_summ_dn2 |> filter(gene_name == "TNFSF4")

# Genes with high KME with each module?
kme_df <- 
    signedKME(vst_dn2, orderMEs(bwnet_dn2$MEs)) |>
    as_tibble(rownames = "gene_id") |>
    pivot_longer(-gene_id, names_to = "module")
    
top_kme_genes <- kme_df |>
    filter(abs(value) > .3) |>
    group_by(module) |>
    top_n(1000, abs(value)) |>
    ungroup() |>
    arrange(module, desc(abs(value)))

top_kme_expression <- top_kme_genes |>
    left_join(select(count_mod_summ_dn2, stim, timep, gene_id, gene_name, vst_counts = value),
	      join_by(gene_id), relationship = "many-to-many") |>
    mutate(module = str_remove(module, "^kME"),
	   g = ntile(module, 2))

mod_dn2_trends_top_p <- 
    top_kme_expression |>
    group_split(g) |>
    map(~ggplot(data = .,
		aes(x = timep, y = vst_counts)) +
	geom_line(aes(group = gene_id), alpha = .1) +
	facet_grid(module~., scale = "free") +
	theme_bw() +
	theme(panel.grid = element_blank()) +
	labs(x = "Hours")) 

ggsave("./plots/wgcna_trends_dn2_top1000.png", 
       cowplot::plot_grid(plotlist = mod_dn2_trends_top_p, ncol = 2),
       height = 7, width = 7)

# GO
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidytext)

run_enrichment <- function(gene_list) {

    enrichGO(gene = gene_list,
	OrgDb = org.Hs.eg.db,
	ont = "BP",
	keyType = "ENSEMBL",
	pAdjustMethod = "fdr",
	pvalueCutoff = 0.05,
	qvalueCutoff = 0.05,
	readable = TRUE) |>
    as.data.frame() |>
    as_tibble()
}

go_res <- 
    top_kme_genes |>
    mutate(gene_id = str_remove(gene_id, "\\.\\d+$"),
	   module = str_remove(module, "^kME")) |>
    {function(x) split(x, x$module)}() |>
    map("gene_id") |>
    map_df(run_enrichment, .id = "module")

go_plot <- go_res |>
    mutate(Description = str_trunc(Description, width = 30)) |>
    group_by(module) |>
    top_n(25, -log10(pvalue)) |>
    ungroup() |>
    ggplot(aes(x = module, 
	       y = reorder_within(Description, by = -log10(pvalue), within = module))) +
	geom_point(aes(fill = -log10(pvalue), size = Count), shape = 21, stroke = .1) +
	scale_y_reordered() +
	scale_size(range = c(0.3, 6)) +
	scale_fill_viridis_c(option = "magma") +
	facet_wrap(~module, nrow = 2, scale = "free") +
	theme_minimal() +
	theme(axis.text.x = element_blank(),
	      axis.text.y = element_text(size = 11),
	      axis.title = element_blank(),
	      panel.grid = element_blank(),
	      legend.position = "top",
	      plot.background = element_rect(fill = "white", color = "white"),
	      strip.text = element_text(size = 14),
	      strip.clip = "off") +
    labs(size = "Number of Genes:", 
	 fill = expression("-log"["10"]("p"))) +
    guides(size = guide_legend(override.aes = list(stroke = 1)),
	   fill = guide_colorbar(barheight = .5, barwidth = 7))

ggsave("./plots/wgcna_dn2_go.png", go_plot, width = 12, height = 10)

# Lupus genes
sle_genes <- "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Khunsriraksakul/Khunsriraksakul_SuppData2.xlsx" |>
    readxl::read_excel(skip = 1) |>
    distinct(`Mapped Gene`) |>
    pull(1)

library(pheatmap)

png("./plots/wgcna_dn2_slegenes.png", units = "in", height = 14, width = 4, res = 200)
tx_to_gene |>
    distinct(gene_id, gene_name) |>
    filter(gene_name %in% sle_genes) |>
    left_join(kme_df) |>
    mutate(module = str_remove(module, "^kME")) |>
    drop_na(module) |> 
    select(-gene_id) |>
    pivot_wider(names_from = module, values_from = value) |>
    column_to_rownames("gene_name") |>
    data.matrix() |>
    pheatmap(fontsize = 9)
dev.off()
