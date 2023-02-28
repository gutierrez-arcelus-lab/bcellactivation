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
    c("grey60", "grey80",
      "grey30", "black",
      brewer.pal(n = 9, "YlOrRd")[c(1, 4)],
      brewer.pal(n = 9, "Blues")[c(2, 8)],
      brewer.pal(n = 9, "Greens")[c(2, 8)],
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
    geom_point(aes(color = demuxlet), size = .25, alpha = .2) +
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
	  plot.margin = margin(t = .5, b = .1, r = .1, l = .1, unit = "cm")) +
    labs(x = "Number of genes detected (in thousands)", 
	 y = "% reads from MT",
	 color = "Genetic demultiplexing:") +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))

ggsave("./plots/qc.png", qcplot, width = 9, height = 2.5)


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
    facet_wrap(~lib, nrow = 1) +
    theme_bw() +
    theme(panel.grid.major.x = element_blank()) +
    labs(x = NULL, y = "Total number of cells")

ggsave("./plots/cells.png", cells_plot, width = 9, height = 2.5)


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

