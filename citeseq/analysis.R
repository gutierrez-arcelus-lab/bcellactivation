# Packages
## single-cell data analysis
library(Seurat)
library(demuxmix)

## Data wrangling
library(dplyr)
library(forcats)
library(purrr)
library(readr)
library(tidyr)

## Plotting
library(ggplot2)
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
    c("Hashtag1" = "BCR 72h",
      "Hashtag2" = "TLR7 72h",
      "Hashtag3" = "BCR 24h", 
      "Hashtag4" = "TLR7 24h",
      "Hashtag5" = "Res 24h",
      "Hashtag6" = "Day 0")

pilot2_stims <- 
    c("Hashtag1" = "Day 0", 
      "Hashtag2" = "IL4 24h",
      "Hashtag3" = "BCR 24h",
      "Hashtag4" = "BCR+TLR7 24h",
      "Hashtag5" = "TLR7 24h", 
      "Hashtag6" = "sCD40L 24h",
      "Hashtag7" = "TLR9 24h",
      "Hashtag8" = "DN2 24h",
      "Hashtag9" = "BCR 72h",
      "Hashtag10" = "BCR+TLR7 72h",
      "Hashtag12" = "TLR7 72h",
      "Hashtag13" = "sCD40L 72h",
      "Hashtag14" = "DN2 72h",
      "Hashtag15" = "TLR9 72h")

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

# Colors
stim_order <- 
    c("Day 0", "Res 24h",
      sprintf("IL4 %sh", c(24, 72)),
      sprintf("sCD40L %sh", c(24, 72)),
      sprintf("BCR %sh", c(24, 72)),
      sprintf("TLR7 %sh", c(24, 72)),
      sprintf("BCR+TLR7 %sh", c(24, 72)),
      sprintf("TLR9 %sh", c(24, 72)),
      sprintf("DN2 %sh", c(24, 72))
    )
      
stim_colors <- 
    c("grey80", "grey60",
      "grey50", "grey40",
      c("goldenrod1", "goldenrod3"),
      brewer.pal(n = 9, "Blues")[c(3, 8)],
      brewer.pal(n = 9, "Greens")[c(3, 8)],
      brewer.pal(n = 9, "Purples")[c(7, 9)],
      grep("pink", colors(), value = TRUE)[c(16, 4)],
      paste0("tomato", c(2, 4))
      )

names(stim_colors) <- stim_order



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

# Genetic demultiplexing
read_demuxlet <- function(f) {
    read_tsv(f) |>
    select(barcode = BARCODE, best = BEST) |>
    extract(best, c("status", "sample"), "([^-]+)-(.+)")
}

demuxlet_pilot2 <- read_demuxlet("./pilot_2/demuxlet/demuxlet_allsnps.best")
demuxlet_1984 <- read_demuxlet("./mgb/demuxlet/demuxlet_results.best") 
demuxlet_1988 <- read_demuxlet("./mgb/demuxlet/demuxlet_1988_results.best")
demuxlet_1990 <- read_demuxlet("./mgb/demuxlet/demuxlet_1990_results.best")

demuxlet_df <- 
    bind_rows("pilot2" = demuxlet_pilot2,
	      "1984" = demuxlet_1984,
	      "1988" = demuxlet_1988,
	      "1990" = demuxlet_1990,
	      .id = "orig.ident") |>
    select(barcode, orig.ident, demuxlet = status)

# Meta data
meta_df <- 
    list(pilot1_obj, pilot2_obj, lib1984_obj, lib1988_obj, lib1990_obj) |>
    map_dfr(function(x) x@meta.data |> 
	    as_tibble(rownames = "barcode") |>
	    select(barcode, orig.ident, n_genes = nFeature_RNA, n_umi = nCount_RNA, percent_mt)) |>
    left_join(demuxlet_df, by = c("barcode", "orig.ident")) |>
    mutate(demuxlet = ifelse(is.na(demuxlet), "UNDEF", demuxlet),
	   orig.ident = fct_inorder(orig.ident))


# QC plots
qcplot <- ggplot(meta_df, aes(x = n_genes, y = percent_mt)) +
    geom_point(aes(color = demuxlet), size = .5, alpha = .5) +
    geom_hline(yintercept = 10, color = "midnightblue", linewidth = 1, alpha = .8) +
    geom_vline(xintercept = 500, color = "midnightblue", linewidth = 1, alpha = .8) +
    scale_x_continuous(labels = function(x) x/1e3 ) +
    scale_color_manual(values = c("UNDEF" = "grey20", "AMB" = "magenta", 
				  "SNG" = "skyblue3", "DBL" = "tomato3")) +
    facet_grid(. ~ orig.ident, scales = "free") +
    theme_bw() +
    theme(panel.grid.minor.x = element_blank(),
	  panel.grid.minor.y = element_blank(),
	  legend.position = "top",
	  plot.margin = margin(t = .25, b = .1, r = .1, l = .1, unit = "cm")) +
    labs(x = "Number of genes detected (in thousands)", 
	 y = "% reads from MT",
	 color = "Genetic demultiplexing:") +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))

# remove bad cells
# compare number of cells

meta_good <- meta_df |>
    filter(n_genes >= 500, percent_mt <= 10) |>
    filter((orig.ident != "pilo1" & demuxlet == "SNG") | (orig.ident == "pilot1"))

good_cells <- split(meta_good, meta_good$orig.ident) |>
    map("barcode")

pilot1_filt <- subset(pilot1_obj, cells = good_cells[["pilot1"]])
pilot2_filt <- subset(pilot2_obj, cells = good_cells[["pilot2"]])
lib1984_filt <- subset(lib1984_obj, cells = good_cells[["1984"]])
lib1988_filt <- subset(lib1988_obj, cells = good_cells[["1988"]])
lib1990_filt <- subset(lib1990_obj, cells = good_cells[["1990"]])

cells_df <- 
    tribble(
	~lib, ~set, ~cells,
	"pilot1", "Before QC", ncol(pilot1_obj),
	"pilot1", "After QC", ncol(pilot1_filt),
	"pilot2", "Before QC", ncol(pilot2_obj),
	"pilot2", "After QC", ncol(pilot2_filt),
	"1984", "Before QC", ncol(lib1984_obj),
	"1984", "After QC", ncol(lib1984_filt),
	"1988", "Before QC", ncol(lib1988_obj),
	"1988", "After QC", ncol(lib1988_filt),
	"1990", "Before QC", ncol(lib1990_obj),
	"1990", "After QC", ncol(lib1990_filt)) |>
    mutate_at(vars(lib, set), fct_inorder)

cells_plot <- ggplot(cells_df, aes(x = set, y = cells)) +
    geom_col(fill = "midnightblue", alpha = .9) +
    scale_y_continuous(labels = function(x) x/1e3) +
    facet_wrap(~lib, nrow = 1) +
    theme_bw() +
    theme(panel.grid.major.x = element_blank()) +
    labs(x = NULL, y = "Total number of cells (k)")


qc_out <- plot_grid(qcplot, NULL, cells_plot, 
		    ncol = 1, rel_heights = c(1, .1, .8),
		    labels = c("A)", "", "B)")) +
    theme(plot.background = element_rect(fill = "white", color = "white"))

ggsave("./plots/qc.png", qc_out, width = 9, height = 5)


# HTO distribution
get_hto <- function(x) {
    x@assays$HTO@counts |> 
	as_tibble(rownames = "hto") |>
	pivot_longer(-hto, names_to = "barcode") |>
	select(barcode, hto, value)
}

hto_counts <- 
    list("pilot1" = pilot1_filt, 
	 "pilot2" = pilot2_filt, 
	 "1984" = lib1984_filt, 
	 "1988" = lib1988_filt,
	 "1990" = lib1990_filt) |>
    map_dfr(get_hto, .id = "orig.ident") |>
    mutate(hto = factor(hto, levels = stim_order),
	   orig.ident = fct_inorder(orig.ident)) |>
    arrange(orig.ident, barcode, hto)

hto_max <- hto_counts |>
    group_by(orig.ident, barcode) |>
    mutate(max_hto = hto[which.max(value)]) |>
    ungroup()

hto_dens <- ggplot(hto_max, aes(x = log10(value + 1))) +
    geom_density(aes(fill = hto), linewidth = .2, alpha = .9) +
    scale_fill_manual(values = stim_colors) +
    facet_wrap(max_hto ~ orig.ident, scales = "free_y", ncol = 6) +
    theme_minimal() +
    theme(text = element_text(size = 10),
	  panel.grid.major.y = element_blank(),
	  panel.grid.minor.y = element_blank(),
	  panel.grid.minor.x = element_blank(),
	  legend.position = "none",
	  strip.text = element_text(size = 8, margin = margin(b = 0.5, t = 0.5)),
	  plot.background = element_rect(fill = "white", color = "white")) +
    labs(x = "Log10 HTO counts",
	 fill = "HTO:")

ggsave("./plots/hto_dens.png", hto_dens, width = 7, height = 7)



# Run HTO demultiplexing
# HTODemux
pilot1_htodemux <- HTODemux(pilot1_filt, assay = "HTO", positive.quantile = 0.99)
pilot2_htodemux <- HTODemux(pilot2_filt, assay = "HTO", positive.quantile = 0.99)
lib1984_htodemux <- HTODemux(lib1984_filt, assay = "HTO", positive.quantile = 0.99)
lib1988_htodemux <- HTODemux(lib1988_filt, assay = "HTO", positive.quantile = 0.99)
lib1990_htodemux <- HTODemux(lib1990_filt, assay = "HTO", positive.quantile = 0.99)

get_htodemux_class <- function(x) {
    x@meta.data |>
	as_tibble(rownames = "barcode") |>
	select(barcode, orig.ident, hto_class = hash.ID, hto_global = HTO_classification.global)
}

htodemux_df <- 
    list(pilot1_htodemux, pilot2_htodemux, lib1984_htodemux, lib1988_htodemux, lib1990_htodemux) |>
    map_dfr(get_htodemux_class)



# demuxmix
run_demuxmix <- function(x) {
    
    hto <- as.matrix(x@assays$HTO@counts)

    rna <- x@meta.data |> 
	select(nFeature_RNA) |>
	as_tibble(rownames = "barcode") |>
	tibble::deframe()

    rna <- rna[colnames(hto)]
    
    dmm <- demuxmix(hto, rna = rna)
    classes <- dmmClassify(dmm)
}


pilot1_dmm <- run_demuxmix(pilot1_filt) |> 
    as_tibble(rownames = "barcode")

pilot2_dmm <- run_demuxmix(pilot2_filt) |> 
    as_tibble(rownames = "barcode")

lib1984_dmm <- run_demuxmix(lib1984_filt) |> 
    as_tibble(rownames = "barcode")

lib1988_dmm <- run_demuxmix(lib1988_filt) |> 
    as_tibble(rownames = "barcode")

lib1990_dmm <- run_demuxmix(lib1990_filt) |> 
    as_tibble(rownames = "barcode")

process_dmm <- function(x) {
    
    x |> 
	mutate(Type = recode(Type, 
			     "multiplet" = "Doublet", 
			     "singlet" = "Singlet", 
			     "negative" = "Negative",
			     "uncertain" = "Uncertain"),
	       HTO = ifelse(Type == "Singlet", HTO, Type)) |>
	select(barcode, hto_class = HTO, hto_global = Type)
}

dmm_df <- 
    list("pilot1" = pilot1_dmm,
	 "pilot2" = pilot2_dmm,
	 "1984" = lib1984_dmm,
	 "1988" = lib1988_dmm,
	 "1990" = lib1990_dmm) |>
    map_dfr(process_dmm, .id = "orig.ident") |>
    select(barcode, orig.ident, everything())


demulti_df <- 
    left_join(htodemux_df, dmm_df, 
	      by = c("barcode", "orig.ident"), 
	      suffix = c(".htodemux", ".dmm")) |>
    mutate_at(vars(hto_global.htodemux, hto_global.dmm), 
	      ~factor(., levels = c("Singlet", "Doublet", "Negative", "Uncertain"))) |>
    mutate(orig.ident = factor(orig.ident, 
			       levels = c("pilot1", "pilot2", "1984", "1988", "1990")))

demulti_summ <- demulti_df |>
    count(orig.ident, hto_global.htodemux, hto_global.dmm) |>
    complete(orig.ident, hto_global.htodemux, hto_global.dmm, fill = list(n = 0))

demux_plot1 <- ggplot(demulti_summ, aes(x = hto_global.htodemux, y = hto_global.dmm)) +
    geom_tile(aes(fill = n), show.legend = FALSE) +
    geom_text(aes(label = n), size = 2.5) +
    scale_fill_gradient(low = "white", high = "slateblue") +
    facet_grid(. ~ orig.ident) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    labs(x = "HTODemux", y = "demuxmix")


demulti_summ2 <- demulti_df |> 
    select(orig.ident, demuxmix = hto_global.dmm, HTODemux = hto_global.htodemux) |>
    pivot_longer(-orig.ident, names_to = "method") |>
    count(orig.ident, method, value) |>
    mutate(method = fct_inorder(method))


demux_plot2 <- ggplot(demulti_summ2, aes(x = value, y = n, fill = method)) +
    geom_col(position = "dodge", alpha = .9) +
    scale_fill_manual(values = c("HTODemux" = "tomato3", "demuxmix" = "midnightblue")) +
    facet_grid(. ~ orig.ident) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
	  legend.position = "bottom",
	  plot.margin = margin(l = 0.66, r = .2, unit = "cm"),
	  panel.grid.major.x = element_blank()) +
    labs(x = NULL, y = "N droplets", fill = NULL)

demux_plot <- 
    plot_grid(demux_plot1, NULL, demux_plot2, 
	      ncol = 1, rel_heights = c(1, .1, 1),
	      labels = c("A)", "", "B)")) +
    theme(plot.background = element_rect(fill = "white", color = "white"))

ggsave("./plots/demux.png", demux_plot, width = 9, height = 5)


# TSNE
DefaultAssay(pilot1_filt) <- "HTO"
DefaultAssay(pilot2_filt) <- "HTO"
DefaultAssay(lib1984_filt) <- "HTO"
DefaultAssay(lib1988_filt) <- "HTO"
DefaultAssay(lib1990_filt) <- "HTO"

pilot1_filt <- pilot1_filt |>
    {function(x) ScaleData(x, features = rownames(x))}() |>
    {function(x) RunPCA(x, features = rownames(x), approx = FALSE)}() |>
    {function(x) RunTSNE(x, dims = 1:nrow(x), perplexity = 100, nthreads = 4)}()

pilot2_filt <- pilot2_filt |>
    {function(x) ScaleData(x, features = rownames(x))}() |>
    {function(x) RunPCA(x, features = rownames(x), approx = FALSE)}() |>
    {function(x) RunTSNE(x, dims = 1:nrow(x), perplexity = 100, nthreads = 4)}()

lib1984_filt <- lib1984_filt |>
    {function(x) ScaleData(x, features = rownames(x))}() |>
    {function(x) RunPCA(x, features = rownames(x), approx = FALSE)}() |>
    {function(x) RunTSNE(x, dims = 1:nrow(x), perplexity = 100, nthreads = 4)}()

lib1988_filt <- lib1988_filt |>
    {function(x) ScaleData(x, features = rownames(x))}() |>
    {function(x) RunPCA(x, features = rownames(x), approx = FALSE)}() |>
    {function(x) RunTSNE(x, dims = 1:nrow(x), perplexity = 100, nthreads = 4)}()

lib1990_filt <- lib1990_filt |>
    {function(x) ScaleData(x, features = rownames(x))}() |>
    {function(x) RunPCA(x, features = rownames(x), approx = FALSE)}() |>
    {function(x) RunTSNE(x, dims = 1:nrow(x), perplexity = 100, nthreads = 4)}()



add_colors <- c("Uncertain" = "turquoise", "Negative" = "beige", "Doublet" = "black")
stim_order_add <- c(names(rev(add_colors)), stim_order)

tsne_df <- 
    list("pilot1" = pilot1_filt, 
	 "pilot2" = pilot2_filt, 
	 "1984" = lib1984_filt, 
	 "1988" = lib1988_filt,
	 "1990" = lib1990_filt) |>
    map_dfr(function(x) x@reductions$tsne@cell.embeddings |>
		as_tibble(rownames = "barcode"),
	    .id = "orig.ident") |>
    mutate(orig.ident = fct_inorder(orig.ident)) |>
    left_join(demulti_df, by = c("barcode", "orig.ident")) |>
    select(barcode, orig.ident, hto_class.htodemux, hto_class.dmm, tSNE_1, tSNE_2) |>
    pivot_longer(hto_class.htodemux:hto_class.dmm, names_to = "method", values_to = "hto_class") |>
    mutate(method = case_when(method == "hto_class.htodemux" ~ "HTODemux",
			      method == "hto_class.dmm" ~ "demuxmix",
			      TRUE ~ NA_character_),
	   method = fct_inorder(method),
	   hto_class = factor(hto_class, levels = stim_order_add))



tsne_plot <- ggplot(tsne_df, aes(x = tSNE_1, y = tSNE_2)) +
    geom_point(aes(color = hto_class), size = .25) +
    scale_color_manual(values = c(stim_colors, add_colors)) +
    facet_grid(orig.ident ~ method) +
    theme_minimal() +
    theme(axis.text = element_blank(),
	  strip.text = element_text(size = 12),
	  panel.grid = element_blank(),
	  plot.background = element_rect(fill = "white", color = "white")) +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    labs(color = NULL, x = "HTO tSNE 1", y = "HTO tSNE 2")

ggsave("./plots/tsne.png", tsne_plot, width = 6.5, height = 10)

# Admixture plots


admix_df <- bind_rows("HTODemux" = htodemux_df, "demuxmix" = dmm_df, .id = "method") |>
    filter(hto_global %in% c("Singlet", "Doublet")) |>
    left_join(hto_counts, by = c("barcode", "orig.ident"), multiple = "all") |>
    group_by(method, barcode, orig.ident) |>
    mutate(top_hto = hto[which.max(value)],
	   p = value/sum(value),
	   p_top = p[which.max(value)]) |>
    ungroup() |>
    mutate(orig.ident = fct_inorder(orig.ident),
	   hto_class = factor(hto_class, levels = stim_order_add),
	   hto_global = factor(hto_global, levels = c("Singlet", "Doublet")))
	

plot_admix <- function(dat) {

    ggplot(dat, aes(x = reorder_within(barcode, by = p_top, within = top_hto), y = p)) +
	geom_col(aes(fill = hto), position = "fill", width = 1.1, show.legend = FALSE) +
	scale_fill_manual(values = stim_colors) +
	scale_y_continuous(expand = c(0, 0)) +
	facet_grid(orig.ident~top_hto, scales = "free_x", space = "free", switch = 'y') +
	theme_minimal() +
	theme(axis.text = element_blank(),
	      strip.text.x = element_blank(),
	      panel.spacing = unit(.1, "lines"),
	      plot.margin = margin(t = 0, b = 0, l = .2, r = .2, unit = "cm"),
	      plot.background = element_rect(fill = "white", color = "white")) +
	labs(x = NULL, y = NULL)
}

## Singlets
# HTOdemux
admix_htodemux_singlet_pilot1 <- admix_df |>
    filter(method == "HTODemux", 
	   orig.ident == "pilot1", 
	   hto_global == "Singlet") |>
    plot_admix()

admix_htodemux_singlet_pilot2 <- admix_df |>
    filter(method == "HTODemux", 
	   orig.ident == "pilot2", 
	   hto_global == "Singlet") |>
    plot_admix()

admix_htodemux_singlet_1984 <- admix_df |>
    filter(method == "HTODemux", 
	   orig.ident == "1984", 
	   hto_global == "Singlet") |>
    plot_admix()

admix_htodemux_singlet_1988 <- admix_df |>
    filter(method == "HTODemux", 
	   orig.ident == "1988", 
	   hto_global == "Singlet") |>
    plot_admix()

admix_htodemux_singlet_1990 <- admix_df |>
    filter(method == "HTODemux", 
	   orig.ident == "1990", 
	   hto_global == "Singlet") |>
    plot_admix()

admix_htodemux_singlets <- 
    plot_grid(NULL,
	      admix_htodemux_singlet_pilot1,
	      admix_htodemux_singlet_pilot2,
	      admix_htodemux_singlet_1984,
	      admix_htodemux_singlet_1988,
	      admix_htodemux_singlet_1990,
	      ncol = 1,
	      rel_heights = c(.25, rep(1, 5)),
	      labels = c("HTODemux", rep("", 5)),
	      label_size = 8) +
    theme(plot.background = element_rect(fill = "white", color = "white"))

# demuxmix
admix_dmm_singlet_pilot1 <- admix_df |>
    filter(method == "demuxmix", 
	   orig.ident == "pilot1", 
	   hto_global == "Singlet") |>
    plot_admix()

admix_dmm_singlet_pilot2 <- admix_df |>
    filter(method == "demuxmix", 
	   orig.ident == "pilot2", 
	   hto_global == "Singlet") |>
    plot_admix()

admix_dmm_singlet_1984 <- admix_df |>
    filter(method == "demuxmix", 
	   orig.ident == "1984", 
	   hto_global == "Singlet") |>
    plot_admix()

admix_dmm_singlet_1988 <- admix_df |>
    filter(method == "demuxmix", 
	   orig.ident == "1988", 
	   hto_global == "Singlet") |>
    plot_admix()

admix_dmm_singlet_1990 <- admix_df |>
    filter(method == "demuxmix", 
	   orig.ident == "1990", 
	   hto_global == "Singlet") |>
    plot_admix()

admix_dmm_singlets <- 
    plot_grid(NULL,
	      admix_dmm_singlet_pilot1,
	      admix_dmm_singlet_pilot2,
	      admix_dmm_singlet_1984,
	      admix_dmm_singlet_1988,
	      admix_dmm_singlet_1990,
	      ncol = 1,
	      rel_heights = c(.25, rep(1, 5)),
	      labels = c("demuxmix", rep("", 5)),
	      label_size = 8) +
    theme(plot.background = element_rect(fill = "white", color = "white"))

ggsave("./plots/admix_singlets.png", 
       plot_grid(admix_htodemux_singlets, admix_dmm_singlets, ncol = 1),
       height = 8, width = 6.5)



## Doublets
# HTOdemux
admix_htodemux_dbl_pilot1 <- admix_df |>
    filter(method == "HTODemux", 
	   orig.ident == "pilot1", 
	   hto_global == "Doublet") |>
    plot_admix()

admix_htodemux_dbl_pilot2 <- admix_df |>
    filter(method == "HTODemux", 
	   orig.ident == "pilot2", 
	   hto_global == "Doublet") |>
    plot_admix()

admix_htodemux_dbl_1984 <- admix_df |>
    filter(method == "HTODemux", 
	   orig.ident == "1984", 
	   hto_global == "Doublet") |>
    plot_admix()

admix_htodemux_dbl_1988 <- admix_df |>
    filter(method == "HTODemux", 
	   orig.ident == "1988", 
	   hto_global == "Doublet") |>
    plot_admix()

admix_htodemux_dbl_1990 <- admix_df |>
    filter(method == "HTODemux", 
	   orig.ident == "1990", 
	   hto_global == "Doublet") |>
    plot_admix()

admix_htodemux_doublets <- 
    plot_grid(NULL,
	      admix_htodemux_dbl_pilot1,
	      admix_htodemux_dbl_pilot2,
	      admix_htodemux_dbl_1984,
	      admix_htodemux_dbl_1988,
	      admix_htodemux_dbl_1990,
	      ncol = 1,
	      rel_heights = c(.25, rep(1, 5)),
	      labels = c("HTODemux", rep("", 5)),
	      label_size = 8) +
    theme(plot.background = element_rect(fill = "white", color = "white"))

# demuxmix
admix_dmm_dbl_pilot1 <- admix_df |>
    filter(method == "demuxmix", 
	   orig.ident == "pilot1", 
	   hto_global == "Doublet") |>
    plot_admix()

admix_dmm_dbl_pilot2 <- admix_df |>
    filter(method == "demuxmix", 
	   orig.ident == "pilot2", 
	   hto_global == "Doublet") |>
    plot_admix()

admix_dmm_dbl_1984 <- admix_df |>
    filter(method == "demuxmix", 
	   orig.ident == "1984", 
	   hto_global == "Doublet") |>
    plot_admix()

admix_dmm_dbl_1988 <- admix_df |>
    filter(method == "demuxmix", 
	   orig.ident == "1988", 
	   hto_global == "Doublet") |>
    plot_admix()

admix_dmm_dbl_1990 <- admix_df |>
    filter(method == "demuxmix", 
	   orig.ident == "1990", 
	   hto_global == "Doublet") |>
    plot_admix()

admix_dmm_doublets <- 
    plot_grid(NULL,
	      admix_dmm_dbl_pilot1,
	      admix_dmm_dbl_pilot2,
	      admix_dmm_dbl_1984,
	      admix_dmm_dbl_1988,
	      admix_dmm_dbl_1990,
	      ncol = 1,
	      rel_heights = c(.25, rep(1, 5)),
	      labels = c("demuxmix", rep("", 5)),
	      label_size = 8) +
    theme(plot.background = element_rect(fill = "white", color = "white"))

ggsave("./plots/admix_doublets.png", 
       plot_grid(admix_htodemux_doublets, admix_dmm_doublets, ncol = 1),
       height = 8, width = 6.5)


