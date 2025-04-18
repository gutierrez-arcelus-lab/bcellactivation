library(tidyverse)
library(AnnotationHub)
library(locuszoomr)
library(rtracklayer)
library(glue)
library(patchwork)

# Figure colors 
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

atac_colors <- stim_colors[c("Unstim 0h", "Unstim 24h", "IL-4c 24h", "TLR7c 24h", "BCRc 24h", "DN2c 24h")]
atac_colors[c("TLR7c 24h", "BCRc 24h", "DN2c 24h")] <- stim_colors[c("TLR7c 24h", "BCRc 48h", "DN2c 48h")] 

# Locus zoom
ah <- AnnotationHub()
ens_data <- ah[["AH98047"]]

gene <- "HLA-B"
window_size <- 0.25e5

loc <- 
    locus(gene = gene,
	  flank = window_size,
	  ens_db = ens_data)

gene_tracks <- 
    gg_genetracks(loc, cex.text = .6) + 
    scale_x_continuous(limits = loc$xrange/1e6,
		       labels = function(x) round(x, 2),
		       expand = c(0, 0)) +
    theme_minimal() +
    theme(plot.background = element_rect(color = "white", fill = "white"))


# ATAC-seq
bigwigs <- 
    "../atacseq/results/bwa/merged_replicate/bigwig" |>
    list.files(full.name = TRUE, pattern = "*.bigWig")

names(bigwigs) <- sub("^([^_]+_\\d+).+$", "\\1", basename(bigwigs))

interv <- GRanges(glue("chr{loc$seqname}"), IRanges(loc$xrange[1], loc$xrange[2]))

atac_ranges <- map(bigwigs, ~import(., which = interv)) 

atac_covered <-
    atac_ranges |>
    map_dfr(~as.data.frame(.) |> as_tibble(), .id = "stim") |>
    dplyr::select(stim, start, end, score)

atac_gaps <- 
    atac_ranges |>
    map_dfr(~ranges(.) |> gaps() |> as.data.frame() |> as_tibble(), .id = "stim") |>
    mutate(score = 0) |>
    dplyr::select(stim, start, end, score)

atac_peaks <- 
    bind_rows(atac_covered, atac_gaps) |>
     mutate(stim = str_replace(stim, "_", " "),
	   stim = str_replace(stim, "unst", "Unstim"),
	   stim = paste0(stim, "h"),
	   stim = str_replace(stim, "IL4", "IL-4c"),
	   stim = str_replace(stim, "TLR7", "TLR7c"),
	   stim = str_replace(stim, "BCR", "BCRc"),
	   stim = str_replace(stim, "DN2", "DN2c"),
	   stim = factor(stim, levels = names(stim_colors))) |>
    arrange(stim, start) |>
    pivot_longer(start:end, names_to = "dummy", values_to = "pos")

atac_plot <- 
    ggplot(atac_peaks) +
    geom_ribbon(aes(x = pos, ymin = 0, ymax = score, color = stim, fill = stim),
		linewidth = .5, outline.type = "full", alpha = .5) +
    scale_x_continuous(limits = loc$xrange, 
		       labels = function(x) round(x/1e6L, 2),
		       expand = c(0, 0)) +
    scale_color_manual(values = atac_colors) +
    scale_fill_manual(values = atac_colors) +
    facet_wrap(~stim, ncol = 1, strip.position = "right") +
    theme_minimal() +
    theme(
	  legend.position = "none",
	  strip.text.y.right = element_text(angle = 0),
	  panel.grid.major.x = element_blank(),
	  panel.grid.minor.x = element_blank(),
	  panel.grid.major.y = element_blank(),
	  panel.grid.minor.y = element_blank(),
	  plot.background = element_rect(color = "white", fill = "white")) +
    labs(x = NULL)
#

plot_out <- 
    atac_plot / gene_tracks + plot_layout(heights = c(1, .3))

ggsave("./test_atac_plot.png", plot_out)
