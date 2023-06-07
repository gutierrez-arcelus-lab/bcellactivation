library(rtracklayer)
library(tidyverse)
library(cowplot)
library(gridExtra)
library(furrr)
library(ggrepel)


plot_susie <- function(region_index) {

    n_gene_rows <- 15
    loc <- regions$locus[region_index]
    res_df <- filter(plot_df, locus == loc)

    tracks_df <- tracks |>
	filter(gwas_locus == unique(res_df$locus)) |>
	select(chr, gene_id, locus = gene_name, feature, i, start, end) |>
	group_by(chr, gene_id, locus) |>
	nest() |>
	ungroup()

    tracks_df$rowid <- rep(1:n_gene_rows, length.out = nrow(tracks_df))
    tracks_df <- unnest(tracks_df, cols = data)

    strand_df <- distinct(tracks, gene_id, locus = gene_name, strand)	

    genes_df <- 
	left_join(tracks_df, strand_df, by = join_by(gene_id, locus)) |>
	group_by(gene_id, locus) |>
	mutate(pos = ifelse(strand == "-", max(end), min(start)) ) |>
	ungroup() |>
	distinct(gene_id, locus, pos) |>
	left_join(distinct(tracks_df, gene_id, locus, rowid), join_by(gene_id, locus))

    arrows_df <- 
	left_join(tracks_df, strand_df, by = join_by(gene_id, locus)) |>
	group_by(gene_id, locus) |>
	mutate(pos = ifelse(strand == "+", max(end), min(start)) ) |>
	ungroup() |>
	distinct(gene_id, locus, strand, x0 = pos) |>
	mutate(x1 = ifelse(strand == "+", x0 + 1e3, x0 - 1e3)) |>
	left_join(select(genes_df, gene_id, locus, rowid), join_by(gene_id, locus))

    max_pip <- res_df |> 
	filter(!is.na(cs), stat == "Susie PIP") |>
	group_by(cs) |>
	slice_max(value) |>
	mutate(stat = paste0(stat, ";GWAS p-value")) |>
	separate_rows(stat, sep = ";")

    pip_plot <- 
	ggplot(data = res_df, 
	       aes(x = pos, y = value)) +
	#geom_vline(data = max_pip, aes(xintercept = pos), linetype = 2, color = "grey35") +
	geom_point() +
	geom_point(data = filter(res_df, !is.na(cs)),
		   aes(color = cs), shape = 21, size = 3) +
	facet_wrap(~fct_rev(stat), ncol = 1, scales = "free_y") +
	theme(axis.text.x = element_blank(),
	      axis.text.y = element_text(margin = margin(0, -.75, 0, 0, unit = "cm")),
	      axis.ticks = element_blank(),
	      panel.grid.minor.y = element_blank(),
	      panel.background = element_rect(fill = "grey96"),
	      legend.position = "top",
	      plot.margin = margin(0, 0.5, 0.5, 0.5, unit = "cm")) +
	labs(x = NULL,
	     y = NULL, 
	     color = "Credible Set:",
	     title = unique(res_df$locus))

    genes_plot <- 
	ggplot(tracks_df |> mutate(rowid = factor(rowid))) +
	#geom_vline(data = max_pip, aes(xintercept = pos), linetype = 2, color = "grey35") +
	geom_segment(data = filter(tracks_df, feature == "intron"),
		     aes(x = start, xend = end, y = rowid, yend = rowid),
		     linewidth = .5, color = "midnightblue") +
	geom_segment(data = filter(tracks_df, feature == "exon"),
		     aes(x = start, xend = end, y = rowid, yend = rowid),
		     linewidth = 3, 
		     color = "midnightblue") +
	geom_segment(data = arrows_df,
		 aes(x = x0, xend = x1, y = rowid, yend = rowid),
		 linewidth = .2,
		 arrow = arrow(length = unit(0.2, "cm")),
		 color = "blue") +
	geom_text_repel(data = genes_df, 
		  aes(x = pos, y = rowid, label = locus),
		  fontface = "italic", 
		  size = 2,
		  nudge_y = .35,
		  direction = "x",
		  min.segment.length = 0,
		  segment.size = .25,
		  force_pull = 100,
		  segment.color = "grey90",
		  color = "black",
		  max.overlaps = 100) +
	scale_x_continuous(limits = range(res_df$pos), labels = function(x) x/1e6L) +
	scale_y_discrete(breaks = 0:n_gene_rows) +
	theme(axis.text = element_blank(),
	      axis.title = element_blank(),
	      axis.ticks = element_blank(),
	      panel.grid = element_blank(),
	      plot.margin = margin(0, 0.5, 0.5, 0.5, unit = "cm"),
	      panel.background = element_rect(fill = "white", color = "white"))

    gr <- 
	map(bigwigs, ~import(., which = interv[region_index])) |>
	map_dfr(~as.data.frame(.) |> as_tibble(), .id = "stim") |>
	mutate(stim = factor(stim, levels = names(stim_colors)))

    gr_df <- gr |>
	mutate(score = score/max(score),
	       bp = map2(start, end, ~.x:.y)) |>
	select(stim, bp, score) |>
	unnest(cols = bp)

    atac <- 
	ggplot(gr_df) +
	geom_line(aes(x = bp, y = score, group = 1, color = stim),
		  linewidth = .4) +
	scale_x_continuous(limits = range(res_df$pos), labels = function(x) x/1e6L) +
	scale_color_manual(values = stim_colors) +
	facet_wrap(~stim, ncol = 1) +
	theme_minimal() +
	theme(axis.ticks.y = element_blank(),
	      axis.text.y = element_blank(),
	      axis.title.y = element_blank(),
	      strip.text = element_blank(),
	      legend.position = "none",
	      panel.grid = element_blank(),
	      plot.margin = margin(0, 0.5, 0.5, 0.5, unit = "cm"),
	      plot.background = element_rect(fill = "white", color = "white")) +
    labs(x = sprintf("Position in %s (Mb)", unique(res_df$chr)))

    p <- 
	plot_grid(get_legend(pip_plot), 
		  pip_plot + theme(legend.position = "none"),
		  genes_plot,
		  atac,
		  ncol = 1, rel_heights = c(.1, 1, .6, .5)) +
	theme(plot.background = element_rect(fill = "white", color = "white"))

    ggsave(sprintf("./plots/susie%s.pdf", loc), p, 
	   width = 7, height = 7
)

}

# Colors
stim_colors <- 
    c("unst_0" = "grey80", 
      "unst_24" = "grey50",
      "IL4_24" = "black", 
      "BCR_24" = "#00a0e4",
      "TLR7_24" = "#488f31",
      "DN2_24" = "#de425b"
      )

## Data
# Genomics windows
regions <- 
    "./data/regions_bentham_1mb.tsv" |>
    read_tsv(col_names = c("region", "locus")) |>
    extract(region, c("chr", "left", "right"), 
	    "(chr[^:]+):(\\d+)-(\\d+)", 
	    convert = TRUE, remove = FALSE)

# GWAS data
summ_stats <- read_tsv("./data/bentham_opengwas_1MbWindows_hg38_summstats.tsv")

# Gene tracks
tracks <- 
    "/lab-share/IM-Gutierrez-e2/Public/vitor/gencode_v38_gene_tracks.tsv" |>
    read_tsv() |>
    group_by(gene_id) |>
    mutate(tss = ifelse(strand == "+", min(start), max(end))) |>
    ungroup() |>
    inner_join(regions, join_by(chr, between(tss, left, right))) |>
    select(chr, gwas_locus = locus, gene_id, gene_name, transcript_id, 
	   strand, feature, i, start, end)

# Analysis
#plan(multisession, workers = availableCores())

pip_df <- read_tsv("./susie_results.tsv") 

plot_df <- 
    left_join(pip_df, summ_stats, 
	      join_by(locus, chr, pos == bp, ref, alt)) |>
    select(locus, chr, pos, logp, pip, cs, coverage) |>
    pivot_longer(logp:pip, names_to = "stat") |>
    mutate(locus = factor(locus, levels = regions$locus),
	   stat = recode(stat, pip = "Susie PIP", logp = "GWAS p-value"),
	   stat = factor(stat, levels = c("Susie PIP", "GWAS p-value")))

# atacseq
interv <- GRanges(regions$chr, IRanges(regions$left, regions$right))

bigwigs <- "../../atacseq/results/bwa/merged_replicate/bigwig" |>
    list.files(full.name = TRUE, pattern = "*.bigWig")

names(bigwigs) <- sub("^([^_]+_\\d+).+$", "\\1", basename(bigwigs))

# Plot
plan(multisession, workers = availableCores())
future_walk(1:nrow(regions), plot_susie)

#ggsave(
#   filename = "./plots/susie.pdf", 
#   plot = marrangeGrob(plot_list, nrow = 1, ncol = 1), 
#   width = 7, height = 7
#)
#
