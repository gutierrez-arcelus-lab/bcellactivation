library(tidyverse)
library(DESeq2)
library(AnnotationHub)
library(locuszoomr)
library(patchwork)
library(rtracklayer)
library(glue)
library(ggforce)
library(ggbeeswarm)
library(cowplot)
library(ggh4x)

select <- dplyr::select
filter <- dplyr::filter
slice <- dplyr::slice

stim_colors <- 
    read_tsv('../figure_colors.txt', col_names = c('stim', 'time', 'color')) |>
    filter(stim %in% c('Unstim', 'IL4', 'TLR7', 'BCR', 'DN2')) |>
    mutate(stim = recode(stim, 
			 'IL4' = 'IL-4c', 
			 'TLR7' = 'TLR7c',
			 'BCR' = 'BCRc', 
			 'DN2' = 'DN2c')) |>
    unite('stim', c(stim, time), sep = ' ') |>
    mutate(stim = paste0(stim, 'h')) |>
    deframe()


# Low-input RNA-seq
dat <- read_rds('../bcell_lowinput/wgcna/data/gene_expression.rds')

sample_table <- 
    tibble(sample_id = colnames(dat$counts)) |>
    separate(sample_id, c('donor_id', 'stim', 'time'), sep = '_', remove = FALSE) |>
    unite('condition', c(stim, time), sep = '_') |>
    column_to_rownames('sample_id')

dds <- 
    DESeqDataSetFromTximport(dat, sample_table, ~1) |>
    estimateSizeFactors() |>
    {function(x) x[, colSums(counts(x)) > 2e6]}() |>
    vst()

tpm_df <- 
    dat$abundance |>
    {function(x) x[grepl('ENSG00000117586\\.\\d+', rownames(x)), rownames(colData(dds)), drop = FALSE]}() |>
    as.data.frame() |>
    as_tibble() |>
    pivot_longer(everything(), names_to = 'sample_id', values_to = 'tpm') |>
    separate(sample_id, c('donor_id', 'stim', 'timep'), sep = '_') |>
    filter(stim %in% c('Unstim', 'IL4', 'TLR7', 'BCR', 'DN2')) |>
    mutate(timep = parse_number(timep)) |>
    mutate(stim = recode(stim, 
			 'IL4' = 'IL-4c', 
			 'TLR7' = 'TLR7c',
			 'BCR' = 'BCRc', 
			 'DN2' = 'DN2c')) |>
    unite('condition', c(stim, timep), sep = ' ', remove = FALSE) |>
    mutate(stim = factor(stim, levels = c('Unstim', 'IL-4c', 'TLR7c', 'BCRc', 'DN2c')),
	   timep = factor(timep, levels = c(0, 4, 24, 48, 72)),
	   condition = paste0(condition, 'h'),
	   condition = factor(condition, levels = names(stim_colors)))

lowinp <- 
    ggplot(tpm_df, aes(x = timep, y = tpm)) +
    geom_quasirandom(aes(fill = condition),
		     method = "smiley", width = .2, 
		     shape = 21, stroke = .1, size = 2) +
    scale_fill_manual(values = stim_colors) +
    facet_grid(.~stim, scale = 'free', space = 'free_x') +
    theme_minimal() +
    theme(axis.text = element_text(size = 8),
	  axis.title = element_text(size = 8),
	  strip.text = element_text(size = 8),
	  panel.spacing.x = unit(0.2, 'lines'),
	  panel.grid.major.y = element_blank(),
	  panel.grid.minor.y = element_blank(),
	  plot.background = element_rect(fill = 'transparent', color = NA)) +
    guides(fill = 'none') +
    labs(x = 'Time (hours)', y = 'TPM')


# Locus zoom + ATAC-seq
gwas_data <-
    "../finemap_and_coloc/finemap/Langefeld/data/summary_stats.tsv" |>
    read_tsv() |>
    filter(grepl("TNFSF4", gene_region)) |>
    select(chr, pos, rsid, alleleA, alleleB, p)

#fine_map <- 
#    "../finemap_and_coloc/finemap/Langefeld/data/compiled_susie_pips.tsv" |>
#    read_tsv() |>
#    filter(grepl("TNFSF4", locus)) |>
#    select(rsid, ref, alt, pip, cs)
#
# Query NCBI to get GRCh38
# use safely and apply function separately to each SNP to avoid breaks
#safe_query <- safely(.f = rsnps::ncbi_snp_query)
#
#ncbi <- 
#    gwas_data |>
#    separate_rows(rsid, sep = ";") |>
#    filter(!is.na(rsid)) |>
#    distinct(rsid) |>
#    pull(rsid) |>
#    map(safe_query)
#
#ncbi_df <- ncbi |> 
#    map_dfr("result") |> 
#    select(rsid = query, chr = chromosome, bp, rsid_new = rsid)
#
#write_rds(ncbi, "./langefeld_ncbi.rds")

ncbi_df <- 
    read_rds("./langefeld_ncbi.rds") |>
    map_dfr("result") |> 
    select(rsid = query, chr = chromosome, bp, rsid_new = rsid)

fine_map <- 
    "../finemap_and_coloc/finemap/Langefeld/data/SuSiE results.SLE EA Immunochip.16OCT2024.xlsx" |>
    readxl::read_excel() |>
    filter(grepl("TNFSF4", region)) |>
    select(rsid = snp, cs, pip)

plot_data <- 
    gwas_data |>
    mutate(chr = as.character(chr)) |>
    rowid_to_column() |>
    separate_rows(rsid, sep = ";") |>
    left_join(ncbi_df, join_by(chr, rsid)) |>
    filter(!is.na(bp)) |>
    group_by(rowid) |>
    slice(1) |>
    ungroup() |>
    select(chr, pos = bp, rsid, alleleA, alleleB, p) |>
    left_join(fine_map) |>
    mutate(pip_interv = case_when(pip > 0 & pip <= 0.2 ~ "0.2",
				  pip > 0.2 & pip <= 0.4 ~ "0.4",
				  pip > 0.4 & pip <= 0.6 ~ "0.6",
				  pip > 0.6 & pip <= 0.8 ~ "0.8",
				  pip > 0.8 & pip <= 1 ~ "1"),
	   pip_interv = factor(pip_interv, levels = c("0.2", "0.4", "0.6", "0.8", "1")))

top_var <- plot_data |> arrange(p) |> slice(1) |> pull(pos)

top_window <- c(top_var - 1e3, top_var + 1e3)

ah <- AnnotationHub()
ens_data <- ah[["AH98047"]]

loc <- 
    locus(data = data.table::as.data.table(plot_data), 
	  gene = 'TNFSF4',
	  flank = c(0.5e5, 1e5),
	  ens_db = ens_data)

gwas_plot <- 
    ggplot(plot_data, aes(x = pos, y = -log10(p))) +
#    geom_vline(xintercept = top_var, 
#	       linetype = 2, 
#	       linewidth = .2, 
#	       color = "black") +
    geom_point(data = ~ filter(., is.na(cs)), 
	       color = "grey60") +
    geom_point(data = ~ filter(., cs == 1), 
	       aes(color = pip), 
	       size = 2) +
    ggrepel::geom_text_repel(data = ~ filter(., cs == 1),
			     aes(label = rsid), 
			     size = 3, box.padding = .5, 
			     nudge_x = -5000, direction = "x",
			     min.segment.length = 0) +
    scale_color_stepsn(name = "PIP:",
		       limits = c(0, 1), 
		       breaks = c(0, .2, .4, .6, .8, 1),
		       colors = c("#FEEBE2", "#FBB4B9", "#F768A1", "#C51B8A", "#7A0177")) +
    scale_x_continuous(limits = loc$xrange,
		       labels = function(x) x/1e6L,
		       expand = c(0, 0)) +
    theme_minimal() +
    labs(tag = "a") +
    theme(panel.grid = element_blank(),
	  axis.title.x = element_blank(),
	  axis.text.x = element_blank(),
	  axis.title.y = element_text(size = 8, margin = margin(r = 0, unit = "lines")),
	  axis.text.y = element_text(size = 8, margin = margin(r = 0, unit = "lines")),
	  plot.tag = element_text(size = 12, face = "bold"),
	  plot.margin = margin(0, 0, 0, 0),
	  legend.text = element_text(size = 8),
	  legend.title = element_text(size = 8),
	  legend.position.inside = c(0.15, 0.7)
	  ) +
    guides(color = guide_colorbar(position = "inside", 
				  direction = "horizontal",
				  barheight = .3,
				  title.position = "top")) +
    coord_cartesian(clip = "off")


## Gene tracks
g <- 
    gg_genetracks(loc, cex.text = .6, 
		  filter_gene_name = loc$TX$symbol[loc$TX$symbol != ""]) +
#    geom_vline(xintercept = top_var/1e6, 
#		 linetype = 2, linewidth = .2, color = "black") +
    scale_x_continuous(limits = loc$xrange/1e6,
		       labels = function(x) x/1e6L,
		       expand = c(0, 0)) +
    theme_minimal() +
    theme(
	  axis.text = element_blank(),
	  axis.title = element_blank(),
	  axis.line.x = element_blank(),
	  axis.ticks.x = element_blank(),
	  panel.grid = element_blank(),
	  plot.margin = margin(0, 0, -1, 0, unit = "lines"))

## ATAC-seq
bigwigs <- 
    "../atacseq/results/bwa/merged_replicate/bigwig" |>
    list.files(full.name = TRUE, pattern = "*.bigWig")

names(bigwigs) <- sub("^([^_]+_\\d+).+$", "\\1", basename(bigwigs))

interv <- 
    GRanges(glue("chr{loc$seqname}"), IRanges(loc$xrange[1], loc$xrange[2]))

gr <- 
    map(bigwigs, ~import(., which = interv)) |>
    map_dfr(~as.data.frame(.) |> as_tibble(), .id = "stim") |>
    mutate(stim = str_replace(stim, "_", " "),
	   stim = str_replace(stim, "unst", "Unstim"),
	   stim = paste0(stim, "h"),
	   stim = str_replace(stim, "IL4", "IL-4c"),
	   stim = str_replace(stim, "TLR7", "TLR7c"),
	   stim = str_replace(stim, "BCR", "BCRc"),
	   stim = str_replace(stim, "DN2", "DN2c"),
	   stim = factor(stim, levels = names(stim_colors)))

gr_df <- 
    gr |>
    mutate(
	   #score = score/max(score),
	   bp = map2(start, end, ~.x:.y)) |>
    select(stim, bp, score) |>
    unnest(cols = bp) 


atac_colors <- stim_colors[c("Unstim 0h", "Unstim 24h", "IL-4c 24h", "TLR7c 24h", "BCRc 24h", "DN2c 24h")]
atac_colors[c("TLR7c 24h", "BCRc 24h", "DN2c 24h")] <- stim_colors[c("TLR7c 24h", "BCRc 48h", "DN2c 48h")] 

atac_plot <- 
    ggplot(gr_df) +
#    geom_vline(xintercept = top_var,
#	       linetype = 2, linewidth = .2, color = "black") +
    geom_ribbon(aes(x = bp, ymin = 0, ymax = score, fill = stim, color = stim), 
		linewidth = .5, outline.type = 'full', alpha = .5) +
    scale_x_continuous(limits = loc$xrange, 
		       labels = function(x) round(x/1e6L, 2),
		       expand = c(0, 0)) +
    scale_color_manual(values = atac_colors) +
    scale_fill_manual(values = atac_colors) +
    facet_wrap(~fct_rev(stim), ncol = 1, strip.position = "right") +
    labs(x = "Chromosome 1 (Mb)", y = "ATAC-seq peak score", tag = "b") +
    theme_minimal() +
    theme(
	  axis.text.x = element_text(size = 8),
	  axis.text.y = element_blank(),
	  axis.title.x = element_text(size = 8),
	  axis.title.y = element_text(size = 8, hjust = .15,
				      margin = margin(r = 0, unit = "lines")),
	  legend.position = "none",
	  strip.text.y.right = element_text(size = 8, angle = 0, 
					    margin = margin(t = 6, l = -6.5, unit = "lines")),
	  strip.clip = "off",
	  panel.grid = element_blank(),
	  panel.spacing.y = unit(-3, "lines"),
	  plot.margin = margin(0, 0, 0, 0, unit = "lines"),
	  plot.tag.position = c(0, .8),
	  plot.tag = element_text(size = 12, face = "bold", hjust = 0, vjust = 0)) +
    annotate(geom = "text", x = 173.14e6, y = 0.75, size = 2.5,
	     label = glue("[0–{round(max(gr_df$score), 2)}]"))

fig_ab_top <- 
    gwas_plot / g / atac_plot + 
    plot_layout(heights = c(.4, .1, 1)) &
    theme(plot.background = element_rect(fill = "transparent", color = NA))

# Fig A and B
fig_ab_title <- 
    ggdraw() + 
    draw_label(
	       "SLE risk variant at the TNFSF4/OX40L locus occurs at a DN2c-specific open chromatin region",
	       x = 0,
	       size = 9,
	       hjust = 0
	       ) +
    theme(text = element_text(size = 9),
	  plot.background = element_rect(fill = "transparent", color = NA),
	  plot.margin = margin(l = 1.25, b = -2.5, unit = "lines"))

fig_ab_grid <- 
    plot_grid(fig_ab_title, 
	      ggdraw() + 
		  draw_line(x = c(0.538, 0.538), y = c(0.1, .95), linetype = 2, linewidth = .25) +
		  draw_plot(fig_ab_top),
	      ncol = 1, rel_heights = c(.1, 1))

#fig_ab_grid <- 
#    plot_grid(fig_ab_title, 
#	      ggdraw(fig_ab_top) +
#		  draw_plot(ggplot(plot_data) + 
#				geom_vline(xintercept = top_var) + 
#				scale_x_continuous(limits = loc$xrange, expand = c(0, 0)) +
#				theme_minimal() +
#				theme(text = element_blank())),				
#	      ncol = 1, rel_heights = c(.1, 1))


# Fig C
fig_c_title <- 
    ggdraw() + 
    draw_label(
	       "TNFSF4/OX40L gene expression is induced in late\ntime points after BCR stimulation",
	       x = 0,
	       size = 9,
	       hjust = 0
	       ) +
    theme(text = element_text(size = 9),
	  plot.margin = margin(l = 1.25, unit = "lines"))

fig_c_grid <- 
    plot_grid(fig_c_title, lowinp, ncol = 1, rel_heights = c(.1, 1))


# Fig D
fig_d_title <- 
    ggdraw() + 
    draw_label(
	       "CRISPR inactivation near rs2205960\naffects TNFSF4 expression",
	       x = 0,
	       size = 9,
	       hjust = 0
	       ) +
    theme(text = element_text(size = 9),
	  plot.margin = margin(l = 1.25, unit = "lines"))


crispr_1 <- 
    readxl::read_excel("./OX40L qPCR raw data.xlsx", 1) |>
    rownames_to_column("sample_id") |>
    pivot_longer(-sample_id, names_to = "guide")

crispr_2 <- 
    readxl::read_excel("./OX40L qPCR raw data.xlsx", 2) |>
    rownames_to_column("sample_id") |>
    pivot_longer(-sample_id, names_to = "guide")
   
crispr_df <- 
    bind_rows("Cas9-KRAB" = crispr_1, "Cas9-KRAB-MeCP2" = crispr_2, .id = "version") |>
    mutate(guide = str_remove(guide, "_")) |>
    mutate_at(vars(version, guide), fct_inorder) |>
    mutate(control_test = ifelse(guide == "SCR", "control", "test"))


fig_d <- 
    ggplot(crispr_df, aes(x = guide, y = value)) +
    geom_quasirandom(aes(fill = control_test),
		     method = "smiley", width = .2, 
		     shape = 21, stroke = 1, size = 2) +
    geom_text(data = tribble(~version, ~x, ~y, ~lab,
			     "Cas9-KRAB", 1.5, 1.25, "****",
			     "Cas9-KRAB", 2, 1.45, "**",
			     "Cas9-KRAB-MeCP2", 1.5, 1.25, "**",
			     "Cas9-KRAB-MeCP2", 2, 1.45, "****"),
	      aes(x = x, y = y, label = lab)) +
    geom_segment(x = 1, xend = 2, y = 1.2, yend = 1.2) +
    geom_segment(x = 1, xend = 3, y = 1.4, yend = 1.4) +
    scale_y_continuous(limits = c(0, 1.5), breaks = c(0, .5, 1)) +
    scale_fill_manual(values = c("control" = "white", "test" = "red")) +
    facet_wrap(~version, nrow = 1) +
    theme_minimal() +
    theme(
	  axis.text.x = element_text(size = 8, angle = 30, hjust = 1, vjust = 1.5),
	  axis.text.y = element_text(size = 8),
	  panel.grid.minor.y = element_blank(),
	  axis.title = element_text(size = 8),
	  strip.text = element_text(size = 8),
	  strip.clip = "off",
	  plot.margin = margin(t = .25, b = -.5, unit = "lines"),
	  plot.background = element_rect(fill = "transparent", color = NA)) +
    labs(x = NULL, 
	 y = expression(paste("  ", italic("TNFSF4"), "/", italic("GAPDH"), " expression"))) +
    guides(fill = "none")

fig_d_grid <- 
    plot_grid(fig_d_title, fig_d, ncol = 1, rel_heights = c(.1, 1))

bottom_panel <- 
    plot_grid(fig_c_grid, NULL, fig_d_grid, nrow = 1, rel_widths = c(1, .05, .66),
	      axis = "t",
	      labels = c("c", "d", ""), label_size = 12, vjust = .5)

ggsave("fig5.png", 
       plot_grid(fig_ab_grid, NULL, bottom_panel, ncol = 1, rel_heights = c(1, .025, .4)) +
	   theme(plot.background = element_rect(fill = "white", color = "white")),
       width = 6.5, height = 6, dpi = 600)

