# Packages
## single-cell data analysis
library(Seurat)
library(sclibr)
library(demuxmix)

## Data wrangling
library(dplyr)
library(forcats)
library(purrr)
library(readr)
library(tidyr)
library(tibble)

## Plotting
library(ggplot2)
library(tidytext)
library(RColorBrewer)
library(scico)
library(ggsci)
library(cowplot)

# RAM
unix::rlimit_as(1e12)


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

pilot2_reseq_dir <- 
    file.path(labshr, "CITEseq_pilot_2",
	      "SN0264956/broad/hptmp/sgurajal/bwh10x/KW10275_mgutierrez/220909_10X_KW10275_bcl",
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

pilot2_reseq_obj <- make_seurat(pilot2_reseq_dir, project_id = "pilot2-reseq", 
				hto_names = pilot2_stims, 
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

demuxlet_pilot2 <- read_demuxlet("./demultiplexing/demuxlet/demuxlet_pilot2.best")
demuxlet_pilot2_reseq <- read_demuxlet("./demultiplexing/demuxlet/demuxlet_pilot2_reseq.best")
demuxlet_1984 <- read_demuxlet("./demultiplexing/demuxlet/demuxlet_1984_results.best") 
demuxlet_1988 <- read_demuxlet("./demultiplexing/demuxlet/demuxlet_1988_results.best")
demuxlet_1990 <- read_demuxlet("./demultiplexing/demuxlet/demuxlet_1990_results.best")

demuxlet_df <- 
    bind_rows("pilot2" = demuxlet_pilot2,
	      "pilot2-reseq" = demuxlet_pilot2_reseq,
	      "1984" = demuxlet_1984,
	      "1988" = demuxlet_1988,
	      "1990" = demuxlet_1990,
	      .id = "orig.ident") |>
    select(barcode, orig.ident, demuxlet = status)

# Meta data
meta_df <- 
    list(pilot1_obj, pilot2_obj, pilot2_reseq_obj, lib1984_obj, lib1988_obj, lib1990_obj) |>
    map_dfr(function(x) x@meta.data |> 
	    as_tibble(rownames = "barcode") |>
	    select(barcode, orig.ident, n_genes = nFeature_RNA, n_umi = nCount_RNA, percent_mt)) |>
    left_join(demuxlet_df, by = c("barcode", "orig.ident")) |>
    mutate(demuxlet = case_when(orig.ident == "pilot1" ~ NA,
				orig.ident != "pilot1" & is.na(demuxlet) ~ "UND",
				orig.ident != "pilot1" & !is.na(demuxlet) ~ demuxlet),
	   orig.ident = fct_inorder(orig.ident))


# QC plots
qcplot <- ggplot(meta_df, aes(x = n_genes, y = percent_mt)) +
    geom_point(aes(color = demuxlet), size = .5, alpha = .5) +
    geom_hline(yintercept = 10, color = "midnightblue", linewidth = 1, alpha = .5) +
    geom_vline(xintercept = 500, color = "midnightblue", linewidth = 1, alpha = .5) +
    scale_x_continuous(labels = function(x) x/1e3 ) +
    scale_color_manual(values = c("UND" = "black", "AMB" = "magenta", 
				  "SNG" = "skyblue3", "DBL" = "tomato3"),
		       na.value = "grey") +
    facet_wrap(~orig.ident, nrow = 2) +
    theme_bw() +
    theme(panel.grid.minor.x = element_blank(),
	  panel.grid.minor.y = element_blank(),
	  legend.position = "top",
	  plot.margin = margin(t = .25, b = .1, r = .1, l = .1, unit = "cm")) +
    labs(x = "Number of genes detected (in thousands)", 
	 y = "% reads from MT",
	 color = "Genetic\ndemux:") +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))

# remove bad cells
# compare number of cells
meta_good <- meta_df |>
    filter(n_genes >= 500, percent_mt <= 10) |>
    filter((orig.ident != "pilot1" & demuxlet == "SNG") | (orig.ident == "pilot1"))

good_cells <- split(meta_good, meta_good$orig.ident) |>
    map("barcode")

pilot1_filt <- subset(pilot1_obj, cells = good_cells[["pilot1"]])
pilot2_filt <- subset(pilot2_obj, cells = good_cells[["pilot2"]])
pilot2_reseq_filt <- subset(pilot2_reseq_obj, cells = good_cells[["pilot2-reseq"]])
lib1984_filt <- subset(lib1984_obj, cells = good_cells[["1984"]])
lib1988_filt <- subset(lib1988_obj, cells = good_cells[["1988"]])
lib1990_filt <- subset(lib1990_obj, cells = good_cells[["1990"]])

cells_df <- 
    tribble(
	~lib, ~set, ~cells,
	"pilot1", "Before", ncol(pilot1_obj),
	"pilot1", "After", ncol(pilot1_filt),
	"pilot2", "Before", ncol(pilot2_obj),
	"pilot2", "After", ncol(pilot2_filt),
	"pilot2-reseq", "Before", ncol(pilot2_reseq_obj),
	"pilot2-reseq", "After", ncol(pilot2_reseq_filt),
	"1984", "Before", ncol(lib1984_obj),
	"1984", "After", ncol(lib1984_filt),
	"1988", "Before", ncol(lib1988_obj),
	"1988", "After", ncol(lib1988_filt),
	"1990", "Before", ncol(lib1990_obj),
	"1990", "After", ncol(lib1990_filt)) |>
    mutate_at(vars(lib, set), fct_inorder)

cells_plot <- ggplot(cells_df, aes(x = set, y = cells)) +
    geom_col(fill = "midnightblue", alpha = .9) +
    scale_y_continuous(labels = function(x) x/1e3) +
    facet_wrap(~lib, nrow = 2) +
    theme_bw() +
    theme(text = element_text(size = 9),
	  panel.grid.major.x = element_blank()) +
    labs(x = "QC", y = "Total number of cells (k)")


qc_out <- plot_grid(plot_grid(get_legend(qcplot), 
			      NULL, 
			      nrow = 1, rel_widths = c(1, .5)),
		    plot_grid(qcplot + theme(legend.position = "none"), 
			      NULL,
			      cells_plot, 
			      nrow = 1, rel_widths = c(1, .05, .6),
			      labels = c("A)", "", "B)")),
		    rel_heights = c(.2, 1), 
		    ncol = 1) +
    theme(plot.background = element_rect(fill = "white", color = "white"))

ggsave("./plots/qc.png", qc_out, width = 7, height = 4)


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
	 "pilot2-reseq" = pilot2_reseq_filt,
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


plot_hto <- function(dat, project, n.rows = 2) {

    dat_project <- filter(dat, orig.ident == project)

    ggplot(dat_project, aes(x = log10(value + 1))) +
	geom_density(aes(fill = hto), linewidth = .2, alpha = .9) +
	scale_y_continuous(expand = c(0, 0), breaks = scales::pretty_breaks(3)) +
	scale_fill_manual(values = stim_colors) +
	facet_wrap(~max_hto, nrow = n.rows) +
	coord_cartesian(ylim = c(0, 2.2)) +
	theme_minimal() +
	theme(text = element_text(size = 10),
	      panel.grid.major.y = element_blank(),
	      panel.grid.minor.y = element_blank(),
	      panel.grid.minor.x = element_blank(),
	      legend.position = "none",
	      strip.text = element_text(size = 8, margin = margin(b = 0.5, t = 0.5)),
	      plot.background = element_rect(fill = "white", color = "white")) +
	labs(x = "Log10 HTO counts",
	     fill = "HTO:", title = project)
}

hto_dens_p1 <- plot_hto(hto_max, "pilot1", 1)
hto_dens_p2 <- plot_hto(hto_max, "pilot2", 3)
hto_dens_p2r <- plot_hto(hto_max, "pilot2-reseq", 3)
hto_dens_1984 <- plot_hto(hto_max, "1984", 2)
hto_dens_1988 <- plot_hto(hto_max, "1988", 2)
hto_dens_1990 <- plot_hto(hto_max, "1990", 2)

hto_dens <- plot_grid(hto_dens_p1, hto_dens_p2, hto_dens_p2r,
		      hto_dens_1984, hto_dens_1988, hto_dens_1990,
		      ncol = 1, rel_heights = c(.6, 1, 1, .75, .75, .75))

ggsave("./plots/hto_dens.png", hto_dens, width = 7, height = 11)


# Total HTO across cells
# Does high throughput lead to more HTO reads?
hto_counts_summ <- hto_counts |>
    filter(orig.ident %in% c("1984", "1988", "1990")) |>
    group_by(orig.ident, barcode) |>
    summarise(value = sum(value)) |>
    ungroup()

hto_summ_plot <- 
    ggplot(hto_counts_summ, aes(x = orig.ident, y = log10(value + 1))) +
	geom_jitter(size = .25, alpha = .1) +
	geom_boxplot(fill = NA, color = "cornflowerblue", linewidth = 1, outlier.color = NA) +
	scale_y_continuous(breaks = 1:5, limits = c(1, 5), expand = c(0, 0)) +
	theme_bw() +
	theme(panel.grid.major.x = element_blank(),
	      panel.grid.minor.y = element_blank()) +
	labs(x = "Sequencing library", y = "log10 HTO raw counts")

ggsave("./plots/hto_summ.png", hto_summ_plot, width = 3, height = 4)



# What about gene expression 
demuxlet_ids <- 
    bind_rows("1984" = demuxlet_1984,
	      "1988" = demuxlet_1988,
	      "1990" = demuxlet_1990,
	      .id = "orig.ident") |>
    inner_join(meta_good, by = c("barcode", "orig.ident", "status" = "demuxlet")) |>
    select(orig.ident, barcode, donor = sample) |>
    mutate(tput = case_when(orig.ident == "1984" ~ "low",
			    orig.ident %in% c("1988", "1990") ~ "high",
			    TRUE ~ NA_character_)) |>
    group_by(tput) |>
    mutate(donor = as.numeric(factor(donor))) |>
    ungroup() |>
    select(orig.ident, barcode, donor)

lib1984_gene <- demuxlet_ids |>
    filter(orig.ident == "1984") |>
    {function(x) split(x, x$donor)}() |>
    map_df(function(x) {
	       lib1984_filt@assays$RNA@data[, x$barcode] |>
	       as_tibble(rownames = "gene_id") |>
	       pivot_longer(-gene_id, names_to = "barcode") |>
	       group_by(gene_id) |>
	       summarise(value = mean(value)) |>
	       ungroup()
	      },
	   .id = "donor")

lib1988_gene <- demuxlet_ids |>
    filter(orig.ident == "1988") |>
    {function(x) split(x, x$donor)}() |>
    map_df(function(x) {
	       lib1988_filt@assays$RNA@data[, x$barcode] |>
	       as_tibble(rownames = "gene_id") |>
	       pivot_longer(-gene_id, names_to = "barcode") |>
	       group_by(gene_id) |>
	       summarise(value = mean(value)) |>
	       ungroup()
	      },
	   .id = "donor")

lib1990_gene <- demuxlet_ids |>
    filter(orig.ident == "1990") |>
    {function(x) split(x, x$donor)}() |>
    map_df(function(x) {
	       lib1990_filt@assays$RNA@data[, x$barcode] |>
	       as_tibble(rownames = "gene_id") |>
	       pivot_longer(-gene_id, names_to = "barcode") |>
	       group_by(gene_id) |>
	       summarise(value = mean(value)) |>
	       ungroup()
	      },
	   .id = "donor")

gene_exp_df <- 
    bind_rows("1984" = lib1984_gene,
	      "1988" = lib1988_gene,
	      "1990" = lib1990_gene,
	      .id = "orig.ident") |>
    arrange(orig.ident, donor) |>
    unite("id", c(orig.ident, donor), sep = "_") |>
    group_by(gene_id) |>
    filter(!all(value == 0)) |>
    ungroup() |>
    pivot_wider(names_from = id, values_from = value) |>
    tibble::column_to_rownames("gene_id")

library(GGally)

pairs_plot <- 
    ggpairs(gene_exp_df,
	    diag  = list(continuous = "blankDiag"),
	    lower = list(continuous = wrap("points", alpha = 0.3, size = 0.25))) +
    theme_bw() +
    theme(panel.grid = element_blank())

ggsave("./plots/pairs.png", pairs_plot) 


gene_exp_df |> 
    as_tibble(rownames = "gene_id") |>
    mutate(d = `1984_1` - `1988_2`) |>
    arrange(desc(abs(d))) |>
    slice(1:5) |>
    left_join(pilot1_features, by = "gene_id")




# Run HTO demultiplexing
# HTODemux
pilot1_htodemux <- HTODemux(pilot1_filt, assay = "HTO", positive.quantile = 0.99)
pilot2_htodemux <- HTODemux(pilot2_filt, assay = "HTO", positive.quantile = 0.99)
pilot2_reseq_htodemux <- HTODemux(pilot2_reseq_filt, assay = "HTO", positive.quantile = 0.99)
lib1984_htodemux <- HTODemux(lib1984_filt, assay = "HTO", positive.quantile = 0.99)
lib1988_htodemux <- HTODemux(lib1988_filt, assay = "HTO", positive.quantile = 0.99)
lib1990_htodemux <- HTODemux(lib1990_filt, assay = "HTO", positive.quantile = 0.99)

get_htodemux_class <- function(x) {
    x@meta.data |>
	as_tibble(rownames = "barcode") |>
	select(barcode, orig.ident, hto_class = hash.ID, hto_global = HTO_classification.global)
}

htodemux_df <- 
    list(pilot1_htodemux, pilot2_htodemux, pilot2_reseq_htodemux, 
	 lib1984_htodemux, lib1988_htodemux, lib1990_htodemux) |>
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

pilot2_reseq_dmm <- run_demuxmix(pilot2_reseq_filt) |> 
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
	 "pilot2-reseq" = pilot2_reseq_dmm,
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
			       levels = c("pilot1", "pilot2", "pilot2-reseq", "1984", "1988", "1990")))

demulti_summ <- demulti_df |>
    count(orig.ident, hto_global.htodemux, hto_global.dmm) |>
    complete(orig.ident, hto_global.htodemux, hto_global.dmm, fill = list(n = 0))

demux_plot1 <- ggplot(demulti_summ, aes(x = hto_global.htodemux, y = hto_global.dmm)) +
    geom_tile(aes(fill = n), show.legend = FALSE) +
    geom_text(aes(label = n), size = 2.5) +
    scale_fill_gradient(low = "white", high = "slateblue") +
    facet_wrap(~orig.ident, nrow = 2) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
	  plot.margin = margin(t = 10)) +
    labs(x = "HTODemux", y = "demuxmix")

demulti_summ2 <- demulti_df |> 
    select(orig.ident, demuxmix = hto_global.dmm, HTODemux = hto_global.htodemux) |>
    pivot_longer(-orig.ident, names_to = "method") |>
    count(orig.ident, method, value) |>
    mutate(method = fct_inorder(method))

demux_plot2 <- ggplot(demulti_summ2, aes(x = value, y = n, fill = method)) +
    geom_col(position = "dodge", alpha = .9) +
    scale_fill_manual(values = c("HTODemux" = "tomato3", "demuxmix" = "midnightblue")) +
    facet_wrap(~orig.ident, nrow = 2) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
	  legend.position = "top",
	  legend.key.width = unit(.25, "cm"),
	  legend.key.height = unit(.25, "cm"),
	  legend.text = element_text(size = 8),
	  legend.background = element_rect(fill = "transparent"),
	  legend.box.background = element_rect(fill = "transparent", color = NA),
	  panel.grid.major.x = element_blank(),
	  legend.margin = margin(0, 0, 0, 0),
	  legend.box.margin = margin(-5, -5, -5, -5)) +
    labs(x = " ", y = "N droplets", fill = NULL)

demux_plot <- 
    plot_grid(demux_plot1, NULL, demux_plot2, 
	      nrow = 1, rel_widths = c(1, .05, .7),
	      labels = c("A)", "", "B)")) +
    theme(plot.background = element_rect(fill = "white", color = "white"))

ggsave("./plots/demux.png", demux_plot, width = 8, height = 3.5)



# TSNE
DefaultAssay(pilot1_filt) <- "HTO"
DefaultAssay(pilot2_filt) <- "HTO"
DefaultAssay(pilot2_reseq_filt) <- "HTO"
DefaultAssay(lib1984_filt) <- "HTO"
DefaultAssay(lib1988_filt) <- "HTO"
DefaultAssay(lib1990_filt) <- "HTO"

run_pca_tsne <- function(x) {
    
    x |>
    {function(x) ScaleData(x, features = rownames(x))}() |>
    {function(x) RunPCA(x, features = rownames(x), approx = FALSE)}() |>
    {function(x) RunTSNE(x, dims = 1:nrow(x), perplexity = 100, nthreads = 4)}()
}

pilot1_filt <- run_pca_tsne(pilot1_filt)
pilot2_filt <- run_pca_tsne(pilot2_filt)
pilot2_reseq_filt <- run_pca_tsne(pilot2_reseq_filt)
lib1984_filt <- run_pca_tsne(lib1984_filt)
lib1988_filt <- run_pca_tsne(lib1988_filt)
lib1990_filt <- run_pca_tsne(lib1990_filt)


add_colors <- c("Uncertain" = "turquoise", "Negative" = "beige", "Doublet" = "black")
stim_order_add <- c(names(rev(add_colors)), stim_order)

tsne_df <- 
    list("pilot1" = pilot1_filt, 
	 "pilot2" = pilot2_filt, 
	 "pilot2-reseq" = pilot2_reseq_filt, 
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

ggsave("./plots/tsne.png", tsne_plot, width = 6.5, height = 11)


# Admixture plots
admix_df <- bind_rows("HTODemux" = htodemux_df, "demuxmix" = dmm_df, .id = "method") |>
    filter(hto_global %in% c("Singlet", "Doublet")) |>
    left_join(hto_counts, by = c("barcode", "orig.ident"), multiple = "all") |>
    mutate(hto = as.character(hto)) |>
    group_by(method, barcode, orig.ident) |>
    mutate(top_hto = hto[which.max(value)],
	   p = value/sum(value),
	   p_top = p[which.max(value)]) |>
    ungroup() |>
    mutate(hto_facet = ifelse(hto_global == "Singlet", hto_class, top_hto),
	   orig.ident = fct_inorder(orig.ident),
	   hto_class = factor(hto_class, levels = stim_order_add),
	   hto_global = factor(hto_global, levels = c("Singlet", "Doublet")),
	   hto_facet = factor(hto_facet, levels = stim_order))
	

plot_admix <- function(dat) {

    ggplot(dat, aes(y = reorder_within(barcode, by = p_top, within = hto_facet), 
		    x = p)) +
	geom_col(aes(fill = hto), position = "fill", width = 1.1, show.legend = FALSE) +
	scale_fill_manual(values = stim_colors) +
	scale_x_continuous(expand = c(0, 0)) +
	facet_grid(hto_facet ~ hto_global, scales = "free", space = "free") +
	theme_minimal() +
	theme(axis.text = element_blank(),
	      strip.text.y = element_blank(),
	      panel.spacing = unit(.1, "lines"),
	      plot.margin = margin(t = .5, b = 0, l = 0.01, r = 0, unit = "cm"),
	      plot.background = element_rect(fill = "white", color = "white")) +
	labs(x = NULL, y = NULL)
}

plot_admix_panel <- function(dat, demuxtool, project) {

    sng <- dat |>
	filter(method == demuxtool, 
	       orig.ident == project, 
	       hto_global == "Singlet")

    sng_plot <- sng |>
	plot_admix() + 
	labs(title = project) + 
	theme(plot.title = element_text(size = 8, hjust = .5)) 

    dbl <- dat |>
	filter(method == demuxtool, 
	       orig.ident == project, 
	       hto_global == "Doublet")

    dbl_plot <- dbl |>
	plot_admix() + 
	labs(title = " ") + 
	theme(plot.title = element_text(size = 8, hjust = .5)) 

	plot_grid(sng_plot, dbl_plot,
		  ncol = 1, rel_heights = c(1, nrow(dbl)/nrow(sng))) +
	theme(plot.background = element_rect(fill = "white", color = "white"))
}

## Singlets
# HTOdemux
admix_htodemux_pilot1 <- plot_admix_panel(admix_df, "HTODemux", "pilot1")
admix_htodemux_pilot2 <- plot_admix_panel(admix_df, "HTODemux", "pilot2")
admix_htodemux_pilot2_reseq <- plot_admix_panel(admix_df, "HTODemux", "pilot2-reseq")
admix_htodemux_1984 <- plot_admix_panel(admix_df, "HTODemux", "1984")
admix_htodemux_1988 <- plot_admix_panel(admix_df, "HTODemux", "1988")
admix_htodemux_1990 <- plot_admix_panel(admix_df, "HTODemux", "1990")

# dmm
admix_dmm_pilot1 <- plot_admix_panel(admix_df, "demuxmix", "pilot1")
admix_dmm_pilot2 <- plot_admix_panel(admix_df, "demuxmix", "pilot2")
admix_dmm_pilot2_reseq <- plot_admix_panel(admix_df, "demuxmix", "pilot2-reseq")
admix_dmm_1984 <- plot_admix_panel(admix_df, "demuxmix", "1984")
admix_dmm_1988 <- plot_admix_panel(admix_df, "demuxmix", "1988")
admix_dmm_1990 <- plot_admix_panel(admix_df, "demuxmix", "1990")

admix_htodemux <- plot_grid(admix_htodemux_pilot1,
			    admix_htodemux_pilot2,
			    admix_htodemux_pilot2_reseq,
			    admix_htodemux_1984,
			    admix_htodemux_1988,
			    admix_htodemux_1990,
			    nrow = 1)

admix_dmm <- plot_grid(admix_dmm_pilot1,
		       admix_dmm_pilot2,
		       admix_dmm_pilot2_reseq,
		       admix_dmm_1984,
		       admix_dmm_1988,
		       admix_dmm_1990,
		       nrow = 1)

admix_out <-
    plot_grid(
	plot_grid(NULL, NULL, NULL, nrow = 1, rel_widths = c(1, 0.1, 1), 
		  labels = c("HTODemux", "", "demuxmix"), hjust = -1.5),
	plot_grid(admix_htodemux, NULL, admix_dmm, nrow = 1, rel_widths = c(1, 0.1, 1)),
	ncol = 1, rel_heights = c(.025, 1)) +
    theme(plot.background = element_rect(fill = "white", color = "white"))
    
ggsave("./plots/admix.png", admix_out, width = 10, height = 12)

# Save data
cells_1984 <- filter(lib1984_dmm, Type == "singlet") |>
    select(barcode, hto = HTO, prob = Prob)

cells_1988 <- filter(lib1988_dmm, Type == "singlet") |>
    select(barcode, hto = HTO, prob = Prob)

cells_1990 <- filter(lib1990_dmm, Type == "singlet") |>
    select(barcode, hto = HTO, prob = Prob)

lib1984_filt <- subset(lib1984_filt, cells = cells_1984$barcode)
lib1988_filt <- subset(lib1988_filt, cells = cells_1988$barcode)
lib1990_filt <- subset(lib1990_filt, cells = cells_1990$barcode)

lib1984_filt@meta.data <- lib1984_filt@meta.data |> 
    as_tibble(rownames = "barcode") |>
    left_join(cells_1984, by = "barcode") |>
    column_to_rownames("barcode")

lib1988_filt@meta.data <- lib1988_filt@meta.data |> 
    as_tibble(rownames = "barcode") |>
    left_join(cells_1988, by = "barcode") |>
    column_to_rownames("barcode")

lib1990_filt@meta.data <- lib1990_filt@meta.data |> 
    as_tibble(rownames = "barcode") |>
    left_join(cells_1990, by = "barcode") |>
    column_to_rownames("barcode")

write_rds(lib1984_filt, "./data/seurat_1984_qced.rds")
write_rds(lib1988_filt, "./data/seurat_1988_qced.rds")
write_rds(lib1990_filt, "./data/seurat_1990_qced.rds")
