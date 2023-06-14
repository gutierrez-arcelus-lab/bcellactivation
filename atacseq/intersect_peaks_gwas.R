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
## Bentham
bentham_regions <- "../colocalization/finemap/data/regions_bentham_1mb.tsv" |>
    read_tsv(col_names = c("region", "locus")) |>
    extract(region, c("chr", "left", "right"), 
	    "(chr[^:]+):(\\d+)-(\\d+)", 
	    convert = TRUE, remove = FALSE)

bentham_stats <- 
    "../colocalization/finemap/data/bentham_opengwas_1MbWindows_hg38_summstats.tsv" |>
    read_tsv()

bentham_lead <- "../colocalization/finemap/data/bentham_leadvars.tsv" |>
    read_tsv()

#regions <- bentham_regions
#gwas_stats <- bentham_stats

## Langefeld
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

bigwigs <- bigwigs[c("unst_24", "DN2_24")]


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

peaks_intersect |> count(locus)


#inner_join(cs_df, da_all, join_by(chr == Chr, between(pos, Start, End))) |>
#    group_by(locus, Geneid) |>
#    filter(! any(contrast %in% c("unst_0vsIL4_24", "unst_0vsunst_24", "unst_24vsIL4_24")) ) |>
#    ungroup() |>
#    filter(locus == unique(locus)[13]) |> print(width = Inf)
#

# plot
plot_susie <- function(loc) {
    
     peaks_df <- peaks_intersect |> 
	filter(locus == loc) |>
	distinct(locus, pos, pip, cs, Start, End)
    
    locus_coords <- 
	c(peaks_df$Start, peaks_df$End) |>
	range() +
	`+`(c(-2e5, 2e5))
    
    n_gene_rows <- 15
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
	summarise(pos = min(bp))  |>
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
	geom_point() +
	geom_point(data = filter(res_df, !is.na(cs)),
		   aes(color = cs), shape = 21, size = 3) +
	facet_wrap(~fct_rev(stat), ncol = 1, scales = "free_y") +
	theme(text = element_text(size = 10),
	      axis.text.x = element_blank(),
	      axis.ticks = element_blank(),
	      panel.grid.minor.y = element_blank(),
	      panel.background = element_rect(fill = "grey96"),
	      legend.position = "top") +
	labs(x = NULL,
	     y = NULL, 
	     color = "Credible Set:",
	     title = unique(res_df$locus))

    genes_plot <- 
	ggplot(tracks_df) +
	geom_line(aes(x = bp, 
		      y = factor(rowid), 
		      group = interaction(gene_id, feature, i),
		      linewidth = feature)) +
	geom_text(data = genes_df, aes(x = pos, y = rowid, label = gene_name),
		  fontface = "italic", hjust = 1, size = 2) +
	scale_linewidth_manual(values = c("exon" = 4, "intron" = .5)) +
	scale_x_continuous(limits = range(res_df$pos), labels = function(x) x/1e6L) +
	scale_y_discrete(breaks = 0:n_gene_rows) +
	theme(axis.text = element_blank(),
	      axis.title = element_blank(),
	      axis.ticks = element_blank(),
	      legend.position = "none",
	      panel.grid = element_blank(),
	      plot.margin = margin(0, 0.5, 0.5, 0.5, unit = "cm"),
	      panel.background = element_rect(fill = "white", color = "white"))

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
		   color = "black", linewidth = .2, linetype = 3) +
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
	      plot.background = element_rect(fill = "white", color = "white")) +
    labs(x = sprintf("Position in %s (Mb)", unique(res_df$chr)))

    p <- pip_plot / genes_plot / atac + plot_layout(heights = c(1, 1, .5))
    
    ggsave(sprintf("./plots/langefeld/susie_%s.pdf", loc), p, 
	   width = 5, height = 6, dpi = 200)
}


walk(unique(peaks_intersect$locus), plot_susie)

