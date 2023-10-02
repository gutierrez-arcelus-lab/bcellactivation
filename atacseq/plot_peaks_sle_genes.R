library(rtracklayer)
library(tidyverse)
library(ggrepel)
library(patchwork)
library(furrr)

# Colors
stim_colors <- 
    c("unst 0" = "grey80", 
      "unst 24" = "grey50",
      "IL4 24" = "black", 
      "BCR 24" = "#00a0e4",
      "TLR7 24" = "#488f31",
      "DN2 24" = "#de425b"
      )

# GWAS
bentham_stats <- 
    "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Bentham/harmonized/26502338-GCST003156-EFO_0002690.h.tsv.gz" |>
    read_tsv(col_types = c("chromosome" = "c")) |>
    select(chr = chromosome, pos = base_pair_location, rsid = variant_id, p_value) |>
    mutate(chr = recode(chr, "23" = "X"))

bentham_lead <- 
    "/lab-share/IM-Gutierrez-e2/Public/vitor/lupus/sle_variants/paper_data/bentham_tab1.tsv" |>
    read_tsv() |>
    select(locus, chr, rsid = snp_id) |>
    left_join(bentham_stats, join_by(chr, rsid)) |>
    select(locus, chr, rsid, chr, pos)

bentham_regions <- 
    bentham_lead |>
    mutate(start = pos - 1e5, end = pos + 1e5) |>
    select(locus, chr, start, end)

bentham_stats_regions <- 
    inner_join(bentham_stats, bentham_regions, 
	       join_by(chr, between(pos, start, end))) |>
    select(locus, chr, pos, rsid, p_value)

# Gene tracks
tracks <- 
    "/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens/gene_tracks_gencodev39.tsv" |>
    read_tsv() |>
    mutate(chr = str_remove(chr, "chr"),
	   gene_id = str_remove(gene_id, "\\.\\d+$")) |>
    group_by(gene_id, gene_name) |>
    nest() |>
    group_by(gene_name) |>
    mutate(idx = seq_len(n())) |>
    ungroup() |>
    mutate(gene_name = ifelse(idx > 1, sprintf("%s (%s)", gene_name, idx), gene_name)) |>
    select(-idx) |>
    unnest(cols = data)

tracks_regions <- 
    inner_join(tracks, select(bentham_regions, locus, chr, begin = start, stop = end),
	       join_by(chr, between(start, begin, stop))) |>
    distinct(locus, gene_id, gene_name) |>
    inner_join(tracks, join_by(gene_id, gene_name))

# Atac-seq peaks
bigwigs <- "../atacseq/results/bwa/merged_replicate/bigwig" |>
    list.files(full.name = TRUE, pattern = "*.bigWig")

names(bigwigs) <- 
    sub("^([^\\.]+).+$", "\\1", basename(bigwigs)) |>
    str_replace("_", " ")

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

da_regions <-
    da_all |>
    mutate(Chr = str_remove(Chr, "chr")) |>
    inner_join(bentham_regions, 
	       join_by(Chr == chr, between(End, start, end))) |>
    select(contrast, locus, Chr, Start, End, pvalue, padj) |>
    filter(grepl("unst_24", contrast)) |>
    group_by(contrast) |>
    nest() |>
    separate(contrast, c("c1", "c2"), sep = "vs") |>
    mutate(stim = ifelse(c1 == "unst_24", c2, c1),
	   stim = str_replace(stim, "_", " "),
	   stim = factor(stim, levels = names(stim_colors))) |>
    select(stim, data) |>
    unnest(cols = data)

# plot
plot_region <- function(loc) {
    
    region_coords <- bentham_regions |>
	filter(locus == loc)

    gwas_df <- bentham_stats_regions |>
	filter(locus == loc)

    gwas_plot <-     
	ggplot(data = gwas_df, 
	       aes(x = pos, y = -log10(p_value))) +
	geom_point(size = 2) +
	scale_x_continuous(labels = function(x) x/1e6L) +
	theme_bw() +
	theme(axis.text.x = element_blank(),
	      axis.ticks.x = element_blank(),
	      panel.grid.major.x = element_blank(),
	      panel.grid.minor.x = element_blank(),
	      panel.grid.major.y = element_line(linewidth = .25),
	      panel.grid.minor.y = element_blank(),
	      panel.border = element_blank(),
	      plot.margin = margin(r = 2, l = 2, unit = "in"),
	      plot.background = element_rect(fill = "white", color = "white")) +
	labs(x = NULL,
	     y = "-log10 (P-value)",
	     title = loc) +
	coord_cartesian(clip = "off")

    tracks_df <- tracks_regions |>
	filter(locus == loc, end > region_coords$start, start < region_coords$end) |>
	mutate(start = ifelse(start < region_coords$start, region_coords$start, start),
	       end = ifelse(end > region_coords$end, region_coords$end, end))

    space_labels <- length(region_coords$start:region_coords$end)/5 
    tracks_range <- 
	tracks_df |>
	group_by(gene_id, gene_name) |>
	summarise(start = min(start) - space_labels,
		  end = max(end)) |>
	ungroup() |>
	mutate(bp = map2(start, end, ~.x:.y)) |>
	select(gene_id, bp) |>
	unnest(bp) |>
	group_by(bp) |>
	filter(n() > 1) |>
	summarise(genes = paste(gene_id, collapse = "/")) |>
	ungroup() |>
	distinct(genes) |>
	rowid_to_column("g") |>
	separate_rows(genes, sep = "/")

    gene_ys <- 
	distinct(tracks_df, gene_id, gene_name) |>
	left_join(tracks_range, join_by(gene_id == genes)) |>
	mutate(g = replace_na(g, 0)) |>
	arrange(g) |>
	group_by(gene_id, gene_name) |>
	slice_max(g) |>
	group_by(g) |>
	mutate(s = seq_len(n())) |>
	ungroup() |>
	mutate(y = g + s) |>
	select(gene_name, y)

    tracks_for_plot <- left_join(tracks_df, gene_ys, join_by(gene_name))

    gene_labels <- tracks_for_plot |>
	group_by(gene_name, y) |>
	summarise(s = min(start)) |>
	ungroup()

    arrow_span <- length(region_coords$start:region_coords$end)/100
    arrows_df <- tracks_for_plot |>
	group_by(gene_name) |>
	mutate(i = ifelse(strand == "+", max(end), min(start))) |>
	ungroup() |>
	distinct(gene_name, start = i, y, strand) |>
	mutate(end = ifelse(strand == "+", start + arrow_span, start - arrow_span))

    gene_plot <- 
	ggplot(tracks_for_plot) +
	geom_segment(aes(x = start, xend = end, y = y, yend = y, linewidth = feature),
		     color = "midnightblue") +
	geom_segment(data = arrows_df,
		     aes(x = start, xend = end, y = y, yend = y),
		     arrow = arrow(length = unit(.2, "cm"), type = "closed"),
		     color = "midnightblue", lineend = "round") +
	geom_text(data = gene_labels, 
		  aes(x = s, y = y, label = gene_name),
		  hjust = 1.25,
		  size = 2, 
		  color = "grey40", 
		  fontface = "italic") +
	scale_x_continuous(limits = c(region_coords$start, region_coords$end),
				    labels = function(x) x/1e6L) +
	scale_y_continuous(limits = c(0, max(gene_ys$y) + 1)) +
	scale_linewidth_manual(values = c("exon" = 2, "intron" = 1)) +
	theme_bw() +
	theme(
	      axis.text.x = element_blank(),
	      axis.ticks.x = element_blank(),
	      axis.text.y = element_blank(),
	      axis.ticks.y = element_blank(),
	      legend.position = "none",
	      panel.grid = element_blank() ,
	      panel.border = element_blank(),
	      plot.background = element_rect(fill = "white", color = "white"),
	      plot.margin = margin(r = 2, l = 2, unit = "in")) +
	labs(x = NULL, y = NULL) +
	coord_cartesian(clip = "off")

    peaks_df <- da_regions |> 
	filter(locus == loc)
    
    interv <- 
	GRanges(paste0("chr", region_coords$chr), 
		IRanges(region_coords$start, region_coords$end))

    gr <- 
	map(bigwigs, ~import(., which = interv)) |>
	map_dfr(~as.data.frame(.) |> as_tibble(), .id = "stim") |>
	mutate(stim = factor(stim, levels = names(stim_colors)))

    gr_df <- gr |>
	mutate(score = score/max(score),
	       bp = map2(start, end, ~.x:.y)) |>
	select(stim, bp, score) |>
	unnest(cols = bp)

    diff_peaks_stars <- 
	inner_join(gr_df, 
	       filter(da_regions, locus == loc),
	       join_by(stim, between(bp, Start, End))) |>
	group_by(stim, Start, End) |>
	slice_max(score) |>
	group_by(stim, Start, End, score) |>
	summarise(bp = median(bp)) |>
	ungroup() |>
	select(stim, bp, score)

    atac <- 
	ggplot(gr_df) +
	geom_line(aes(x = bp, y = score, group = 1, color = stim),
		  linewidth = .7) +
	scale_x_continuous(limits = c(region_coords$start, region_coords$end), 
			   labels = function(x) round(x/1e6L, 2)) +
	scale_y_continuous(limits = c(0, 1.1)) +
	scale_color_manual(values = stim_colors) +
	geom_text(data = diff_peaks_stars,
		  aes(x = bp, y = score + .05),
		  size = 4, label = "*", fontface = "bold") +
	facet_grid(stim~., switch = "y") +
	theme_bw() +
	theme(
	      axis.text.y = element_blank(),
	      axis.ticks.y = element_blank(),
	      strip.text.y.left = element_text(angle = 0, ),
	      strip.background = element_rect(fill = "white", color = "white"),
	      strip.placement = "outside",
	      legend.position = "none",
	      panel.grid = element_blank(),
	      panel.border = element_blank(),
	      plot.margin = margin(r = 2, l = 2, unit = "in"),
	      plot.background = element_rect(fill = "white", color = "white")) +
    labs(x = sprintf("Position on chr%s (Mb, GRCh38)", unique(gwas_df$chr)),
	 y = NULL)
    
    ggsave(sprintf("./plots/atac_gwas/bentham_%s.png", loc),
	   plot_spacer() + gwas_plot + gene_plot + atac + plot_spacer() + plot_layout(heights = c(.8, 1, 1, 1, .8), ncol = 1),
	   height = 8.5, width = 11, dpi = 300)

}


walk(bentham_regions$locus, plot_region)




# Try Gviz
library(Gviz)

gtf <- 
    "/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens/gencode.v42.primary_assembly.annotation.gtf" |>
    read_tsv()



gene_models <- 
    tracks_df |>
    filter(feature == "exon") |>
    mutate(width = map2_int(start, end, ~length(.x:.y) - 1L)) |>
    select(chromosome = chr, start, end, strand, gene = gene_id, exon = index, 
	   transcript = transcript_id, symbol = gene_name)

grtrack <- 
    GeneRegionTrack(gene_models, 
		    genome = "hg38", 
		    chromosome = unique(gene_models$chromosome))


displayPars(grtrack) <- 
    list(fill = "midnightblue", col.line = "midnightblue", col = "midnightblue", 
	 min.height = 50, arrowFeatherWidth = 15,  
	 showTitle = FALSE)

png("./plots/test.png", res = 600, width = 7, height = 2, unit = "in")
plotTracks(grtrack, transcriptAnnotation = "symbol")
dev.off()
