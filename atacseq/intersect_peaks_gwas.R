library(rtracklayer)
library(tidyverse)
library(ggrepel)
library(patchwork)
library(furrr)

# Colors
stim_colors <- 
    c("unst_0" = "grey80", 
      "unst_24" = "grey50",
      "IL4_24" = "black", 
      "BCR_24" = "#00a0e4",
      "TLR7_24" = "#488f31",
      "DN2_24" = "#de425b"
      )

# GWAS data
# Bentham
#bentham_regions <- "../colocalization/finemap/data/regions_bentham_1mb.tsv" |>
#    read_tsv(col_names = c("region", "locus")) |>
#    extract(region, c("chr", "left", "right"), 
#	    "(chr[^:]+):(\\d+)-(\\d+)", 
#	    convert = TRUE, remove = FALSE)
#
#bentham_stats <- 
#    "../colocalization/finemap/data/bentham_opengwas_1MbWindows_hg38_summstats.tsv" |>
#    read_tsv()
#
#bentham_lead <- "../colocalization/finemap/data/bentham_leadvars.tsv" |>
#    read_tsv()
#
#regions <- bentham_regions
#gwas_stats <- bentham_stats

# Langefeld
langefeld_regions <- "../colocalization/finemap/data/regions_langefeld_1mb.tsv" |>
    read_tsv(col_names = c("region", "locus")) |>
    extract(region, c("chr", "left", "right"), 
	    "(chr[^:]+):(\\d+)-(\\d+)", 
	    convert = TRUE, remove = FALSE)

langefeld_stats <- 
    "../colocalization/finemap/data/langefeld_1MbWindows_hg38_summstats.tsv" |>
    read_tsv() |>
    mutate(logp = -log10(p)) |>
    select(chr, locus, varid, bp, ref = major, alt = minor, beta, se, logp)

langefeld_lead <- 
    "../sle_variants/paper_data/langefeld_top.tsv" |>
    read_tsv()

regions <- langefeld_regions
gwas_stats <- langefeld_stats

# Gene tracks
tracks <- 
    "/lab-share/IM-Gutierrez-e2/Public/vitor/genecode_V38_tracks.tsv" |>
    read_tsv() |>
    group_by(gene_id) |>
    mutate(tss = ifelse(strand == "+", min(start), max(end))) |>
    ungroup() |>
    inner_join(regions, join_by(chr, between(tss, left, right))) |>
    select(chr, gwas_locus = locus, gene_id, gene_name, transcript_id, 
	   strand, feature, i, start, end)

# Susie
#pip_df <- "../colocalization/finemap/susie_results.tsv" |>
#    read_tsv() 
#
pip_df <- "../colocalization/finemap/susie_results_langefeld.tsv" |>
    read_tsv()

cs_df <- filter(pip_df, !is.na(cs))

plot_df <- 
    left_join(pip_df, gwas_stats, 
	      join_by(locus, chr, pos == bp, ref, alt)) |>
    select(locus, chr, pos, logp, pip, cs, coverage) |>
    pivot_longer(logp:pip, names_to = "stat") |>
    mutate(locus = factor(locus, levels = regions$locus),
	   stat = recode(stat, pip = "Susie PIP", logp = "GWAS p-value"),
	   stat = factor(stat, levels = c("Susie PIP", "GWAS p-value")))

# Atac-seq peaks
bigwigs <- "../atacseq/results/bwa/merged_replicate/bigwig" |>
    list.files(full.name = TRUE, pattern = "*.bigWig")

names(bigwigs) <- sub("^([^_]+_\\d+).+$", "\\1", basename(bigwigs))

#bigwigs <- bigwigs[c("unst_24", "DN2_24")]


# Differential accessibility analysis
da_files <- 
    list.files("./results_deseq2", 
	       pattern = "\\.deseq2\\.FDR0\\.05\\.results\\.txt",
	       full.names = TRUE)

da_names <- 
    basename(da_files) |>
    strsplit("\\.") |>
    map_chr(1)

da_all <- da_files |>
    setNames(da_names) |>
    {function(x) x[!grepl("72", x)]}() |>
    map_dfr(~read_tsv(.) |> select(Geneid:padj), .id = "contrast")

peaks_intersect <- 
    inner_join(cs_df, filter(da_all, contrast == "unst_24vsDN2_24"), 
	       join_by(chr == Chr, between(pos, Start, End))) |>
    select(-contrast, -Geneid)



# plot
plot_susie <- function(loc) {
    
     peaks_df <- peaks_intersect |> 
	filter(locus == loc) |>
	distinct(locus, pos, pip, cs, Start, End)
    
    locus_coords <- 
	c(peaks_df$Start, peaks_df$End) |>
	range() +
	`+`(c(-5e4, 5e4))
    
    n_gene_rows <- 5
    res_df <- filter(plot_df, locus == loc, between(pos, locus_coords[1], locus_coords[2]))
    
    interv <- GRanges(unique(res_df$chr), IRanges(min(res_df$pos), max(res_df$pos)))

    tracks_df <- tracks |>
	filter(gwas_locus == unique(res_df$locus)) |>
	mutate(bp = map2(start, end, ~.x:.y)) |>
	select(chr, gene_id, gene_name, feature, strand, i, bp,) |>
	unnest(cols = bp) |>
	filter(between(bp, min(res_df$pos), max(res_df$pos))) |>
	group_by(gene_id) |>
	nest() |>
	ungroup() |>
	mutate(rowid = rep(1:n_gene_rows, length.out = n())) |>
	unnest(cols = data)

    genes_df <- 
	tracks_df |>
	group_by(gene_id, gene_name, rowid) |>
	summarise(pos = median(bp))  |>
	ungroup()

    max_pip <- res_df |> 
	filter(!is.na(cs), stat == "Susie PIP") |>
	group_by(cs) |>
	slice_max(value) |>
	mutate(stat = paste0(stat, ";GWAS p-value")) |>
	separate_rows(stat, sep = ";")

    pip_plot <- 
	ggplot(data = res_df, 
	       aes(x = pos, y = value)) +
	geom_point(size = .5) +
	geom_point(data = filter(res_df, !is.na(cs)),
		   aes(color = cs), shape = 21, size = 1.5) +
	facet_wrap(~fct_rev(stat), ncol = 1, scales = "free_y") +
	theme_minimal() +
	theme(axis.text.x = element_blank(),
	      axis.text.y = element_text(size = 6),
	      axis.ticks = element_blank(),
	      panel.grid = element_blank(),
	      panel.background = element_rect(fill = "white", color = "white"),
	      strip.text = element_text(size = 8),
	      legend.position = "top",
	      legend.margin = margin(0,0,0,0),
	      legend.text = element_text(size = 6),
	      legend.title = element_text(size = 8),
	      legend.spacing.x = unit(0.05, 'cm'),
	      plot.title = element_text(size = 8, margin=margin(t=0, b=-2.5))) +
	labs(x = NULL,
	     y = NULL, 
	     color = "CS:",
	     title = unique(res_df$locus)) +
	coord_cartesian(clip = "off")

    genes_plot <- 
	ggplot(tracks_df) +
	geom_line(aes(x = bp, 
		      y = factor(rowid), 
		      group = interaction(gene_id, feature, i),
		      linewidth = feature)) +
	geom_text(data = genes_df, aes(x = pos, y = rowid, label = gene_name),
		  fontface = "italic", hjust = 1, size = 1.75, nudge_y = 0.5) +
	scale_linewidth_manual(values = c("exon" = 2, "intron" = .5)) +
	scale_x_continuous(limits = range(res_df$pos), labels = function(x) x/1e6L) +
	scale_y_discrete(breaks = 0:n_gene_rows) +
	theme(axis.text = element_blank(),
	      axis.title = element_blank(),
	      axis.ticks = element_blank(),
	      legend.position = "none",
	      panel.grid = element_blank(),
	      panel.background = element_rect(fill = "white", color = "white")) +
	coord_cartesian(clip = "off")

    gr <- 
	map(bigwigs, ~import(., which = interv)) |>
	map_dfr(~as.data.frame(.) |> as_tibble(), .id = "stim") |>
	mutate(stim = factor(stim, levels = names(stim_colors)))

    gr_df <- gr |>
	mutate(score = score/max(score),
	       bp = map2(start, end, ~.x:.y)) |>
	select(stim, bp, score) |>
	unnest(cols = bp)

    atac <- 
	ggplot(gr_df) +
	geom_vline(data = peaks_df, aes(xintercept = pos), 
		   color = "black", linewidth = .25, linetype = 3, alpha = .8) +
	geom_line(aes(x = bp, y = score, group = 1, color = stim),
		  linewidth = .7) +
	scale_x_continuous(limits = range(res_df$pos), 
			   labels = function(x) round(x/1e6L, 2),
			   breaks = scales::pretty_breaks(3)) +
	scale_color_manual(values = stim_colors) +
	facet_wrap(~stim, ncol = 1) +
	theme_minimal() +
	theme(axis.ticks.y = element_blank(),
	      axis.text.x = element_text(size = 6),
	      axis.text.y = element_blank(),
	      axis.title.x = element_text(size = 8),
	      axis.title.y = element_blank(),
	      strip.text = element_blank(),
	      legend.position = "none",
	      panel.grid = element_blank(),
	      plot.background = element_rect(fill = "white", color = "white")) +
    labs(x = sprintf("Position in %s (Mb)", unique(res_df$chr)))

    p <- pip_plot / genes_plot / atac + plot_layout(heights = c(1, .35, .5))
   
    ggsave(sprintf("./plots/langefeld/susie_%s.pdf", loc), p, 
	   width = 5, height = 6, dpi = 200)
}


walk(unique(peaks_intersect$locus), plot_susie)

loc = "ITGAM"
loc = "IL10"
loc = "TNFSF4-LOC100506023"
#
itgam <- p
il10 <- p
tnf <- p

ggsave("./plots/candidates.png",
       cowplot::plot_grid(itgam, tnf, il10, nrow = 1),
       width = 5, height = 3.5)

# Poster
#loc <- "IL10"
loc <- "STAT4"
#loc <- "TNFSF4"

top <- plot_df |>
    filter(grepl(loc, locus), stat == "GWAS p-value") |>
    filter(value == max(value)) |>
    select(locus, chr, pos) |>
    mutate(start = pos - 1e5, end = pos + 4e4)

n_gene_rows <- 5
res_df <- plot_df |>
    filter(grepl(loc, locus), between(pos, top$start, top$end),
	   stat == "GWAS p-value") |>
    select(-stat)

interv <- GRanges(unique(res_df$chr), IRanges(min(res_df$pos), max(res_df$pos)))

diff_access_df <- da_all |>
    inner_join(top, join_by(Chr == chr, between(Start, start, end))) |>
    select(contrast, locus, Chr, Start, End, log2FoldChange, pvalue, padj) |>
    filter(grepl("unst_24", contrast)) |>
    separate(contrast, c("c0", "c1"), sep = "vs") |>
    mutate(stim = ifelse(c0 == "unst_24", c1, c0),
	   stim = factor(stim, levels = names(stim_colors))) |>
    select(stim, locus, Start, End)


tracks_df <- tracks |>
    filter(gwas_locus == unique(res_df$locus)) |>
    mutate(bp = map2(start, end, ~.x:.y)) |>
    select(chr, gene_id, gene_name, feature, strand, i, bp,) |>
    unnest(cols = bp) |>
    filter(between(bp, top$start, top$end)) |>
    group_by(gene_id) |>
    nest() |>
    ungroup() |>
    mutate(rowid = rep(1:n_gene_rows, length.out = n())) |>
    unnest(cols = data)

genes_df <- 
    tracks_df |>
    group_by(gene_id, gene_name, rowid) |>
    summarise(pos = median(bp))  |>
    ungroup()

gwas_plot <- 
    ggplot(data = res_df, 
	   aes(x = pos, y = value)) +
    geom_point(size = 1) +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
	  axis.text.y = element_text(size = 6),
	  axis.title.y = element_text(size = 7),
	  axis.ticks = element_blank(),
	  panel.grid = element_blank(),
	  panel.background = element_rect(fill = "white", color = "white")) +
    labs(x = NULL,
	 y = "GWAS log10 p") +
    coord_cartesian(clip = "off")

genes_plot <- 
    ggplot(tracks_df) +
    geom_line(aes(x = bp, 
		  y = factor(rowid), 
		  group = interaction(gene_id, feature, i),
		  linewidth = feature)) +
    geom_text(data = genes_df, aes(x = pos, y = rowid, label = gene_name),
	      fontface = "italic", hjust = 1, size = 2, nudge_y = 0.5) +
    scale_linewidth_manual(values = c("exon" = 2, "intron" = .5)) +
    scale_x_continuous(limits = range(res_df$pos), labels = function(x) x/1e6L) +
    scale_y_discrete(breaks = 0:n_gene_rows, expand = c(0, 0)) +
    theme(axis.text = element_blank(),
	  axis.title = element_blank(),
	  axis.ticks = element_blank(),
	  legend.position = "none",
	  panel.grid = element_blank(),
	  panel.background = element_rect(fill = "white", color = "white")) +
    coord_cartesian(clip = "off")

gr <- 
    map(bigwigs, ~import(., which = interv)) |>
    map_dfr(~as.data.frame(.) |> as_tibble(), .id = "stim") |>
    mutate(stim = factor(stim, levels = names(stim_colors)))

gr_df <- gr |>
    mutate(score = score/max(score),
	   bp = map2(start, end, ~.x:.y)) |>
    select(stim, bp, score) |>
    unnest(cols = bp)

da_peaks_df <- 
    inner_join(gr_df, diff_access_df, 
	       join_by(stim, between(bp, Start, End))) |>
    group_by(locus, stim, Start, End) |>
    slice_max(score) |>
    group_by(locus, stim, Start, End, score) |>
    summarise(bp = median(bp)) |>
    ungroup()

atac <- 
    ggplot(gr_df) +
    geom_text(data = da_peaks_df, 
	      aes(x = bp, y = score),
	      size = 2.5, label = "*") +
    geom_line(aes(x = bp, y = score, group = 1, color = stim),
	      linewidth = .5) +
    geom_text(data = gr_df |> mutate(bp = min(bp), score = .5) |> distinct(),
	      aes(x = bp, y = score, label = stim, color = stim),
	      size = 2, hjust = "outward") +
    scale_x_continuous(limits = range(res_df$pos), 
		       labels = function(x) round(x/1e6L, 3),
		       breaks = scales::pretty_breaks(3)) +
    scale_color_manual(values = stim_colors) +
    facet_wrap(~stim, ncol = 1) +
    theme_minimal() +
    theme(axis.ticks.y = element_blank(),
	  axis.text.x = element_text(size = 6),
	  axis.text.y = element_blank(),
	  axis.title.x = element_text(size = 8),
	  axis.title.y = element_blank(),
	  strip.text = element_blank(),
	  legend.position = "none",
	  panel.grid = element_blank(),
	  panel.spacing = unit(.05, "cm"),
	  plot.background = element_rect(fill = "white", color = "white")) +
    labs(x = sprintf("Position in %s (Mb)", unique(res_df$chr))) +
    coord_cartesian(clip = "off")

p <- gwas_plot / genes_plot / atac + plot_layout(heights = c(.4, .25, 1))

#il10 <- p
stat4 <- p
#tnf <- p


ggsave("./plots/langefeld/susie_poster_atac.png", 
       cowplot::plot_grid(tnf, il10, stat4, nrow = 1), 
       width = 8, height = 4, dpi = 600)

