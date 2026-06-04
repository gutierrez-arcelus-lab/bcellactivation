library(tidyverse)
library(patchwork)
library(AnnotationHub)
library(ensembldb)
library(locuszoomr)

ah <- AnnotationHub()
ens_data <- ah[["AH98047"]]

conflicted::conflicts_prefer(dplyr::select)
conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::slice)

bigwigs <- 
    list.files("../../shinybcells/inst/extdata/",
	       pattern = "*_filtered.bigWig",
	       full.names = TRUE)
        
names(bigwigs) <- sub("^([^_]+_\\d+).+$", "\\1", basename(bigwigs))

load("../../shinybcells/R/sysdata.rda")

da_files <- 
    list.files("../03_atacseq/2-differential_peaks/results",
               pattern = ".+vs.+\\.tsv$",
               full.names = TRUE) |>
    {function(x) keep(x, grepl("unst_24", x))}()

da_files <- setNames(da_files, basename(da_files) |> str_remove("\\.tsv"))

da_data <-
    map_dfr(da_files, read_tsv, .id = "comparison") |>
    filter(!is.na(padj)) |>
    separate(comparison, c("stim1", "stim2"), sep = "vs")

da_data_fdr <- 
    da_data |>
    filter(padj <= 0.01) |>
    mutate(stim = case_when(stim1 != "unst_24" & stim2 != "unst_0" ~ stim1,
                            TRUE ~ stim2),
           middle = round((Start + End)/2)) |>
    select(stim, interval, Chr, middle)

make_atac_plot <- function(i) {

    loc <-
    	locus(gene = input |> slice(i) |> pull(atac_gene),
    	      flank = input |> slice(i) |> pull(atac_window),
    	      ens_db = ens_data)
    
    interv <-
    	GRanges(paste0("chr", loc$seqname),
    	        IRanges(loc$xrange[1], loc$xrange[2]))

    atac_ranges <- 
    	bigwigs |>
    	map(~rtracklayer::import(., which = interv))

    atac_covered <-
    	atac_ranges |>
    	map_dfr(as.data.frame, .id = "stim") |>
    	select(stim, start, end, score)
    
    atac_gaps <-
    	atac_ranges |>
    	map_dfr(~keepSeqlevels(., paste0("chr", loc$seqname)) |> 
    	            gaps(start = loc$xrange[1], end = loc$xrange[2]) |> 
    	            as.data.frame(),
    	        .id = "stim"
    	) |>
    	mutate(score = 0) |>
    	select(stim, start, end, score)
			
    atac_peaks <-
    	bind_rows(atac_covered, atac_gaps) |>
    	mutate(stim = str_replace(stim, "unst", "Unstim"),
    	       stim = str_replace(stim, "IL4", "IL-4c"),
    	       stim = str_replace(stim, "TLR7", "TLR7c"),
    	       stim = str_replace(stim, "BCR", "BCRc"),
    	       stim = str_replace(stim, "DN2", "DN2c"),
    	       stim = factor(stim, levels = names(atac_colors))) |>
    	arrange(stim, start) |>
    	pivot_longer(start:end, names_to = "dummy", values_to = "pos")
    
    atac_da <-
        da_data_fdr |>
        filter(Chr == paste0("chr", loc$seqname)) |>
        inner_join(atac_covered, join_by(stim, between(middle, start, end))) |>
    	mutate(stim = str_replace(stim, "unst", "Unstim"),
    	       stim = str_replace(stim, "IL4", "IL-4c"),
    	       stim = str_replace(stim, "TLR7", "TLR7c"),
    	       stim = str_replace(stim, "BCR", "BCRc"),
    	       stim = str_replace(stim, "DN2", "DN2c"),
    	       stim = factor(stim, levels = names(atac_colors)))
    
    plot_atac <-
    	ggplot(atac_peaks) +
    	geom_ribbon(aes(x = pos, ymin = 0, ymax = score, color = stim, fill = stim),
    	            linewidth = .25, outline.type = "full", alpha = .5) +
        geom_text(data = atac_da,
                  aes(x = middle, y = score, label = "*"),
                  nudge_y = 0.5, size = 10, size.unit = "pt") +
    	scale_x_continuous(limits = loc$xrange,
    	                   labels = function(x) round(x/1e6L, 2),
    	                   expand = c(0, 0)) +
    	scale_y_continuous(breaks = round(c(0, max(atac_peaks$score)), 1)) +
    	scale_color_manual(values = atac_colors) +
    	scale_fill_manual(values = atac_colors) +
    	facet_wrap(~stim, ncol = 1, strip.position = "right") +
    	theme_minimal() +
    	theme(
    	      axis.text.x = element_blank(),
    	      axis.text.y = element_text(size = 6, vjust = c(.5, 1.5)),
    	      legend.position = "none",
    	      strip.text.y.right = element_text(angle = 0, size = 6),
    	      panel.grid = element_blank(),
    	      panel.spacing = unit(0, "lines"),
    	      plot.margin = margin(b = 0),
    	      plot.background = element_rect(color = "white", fill = "white")
    	) +
    	labs(x = NULL)
    
    gene_tracks <-
    	gg_genetracks(loc, cex.text = .5) +
    	scale_x_continuous(limits = loc$xrange/1e6,
    			   labels = function(x) round(x, 2),
    			   expand = c(0, 0)
    	) +
    	theme_minimal() +
    	theme(
    	      panel.grid = element_blank(),
    	      axis.text = element_text(size = 6),
    	      axis.title = element_text(size = 6),
    	      plot.margin = margin(t = 0),
    	      plot.background = element_rect(color = "white", fill = "white")
    	)

    plot_atac / gene_tracks + plot_layout(heights = c(1, .2))
}

input <- 
    tribble(~atac_gene, ~atac_window,
	    "CD19", 1.5e4,
	    "CD86", 3e4,
	    "BATF", 3e4,
	    "TBX21", 2.5e4,
	    "FCER2", 1e4,
	    "IGHG4", 1e4
	    )

atac_plot_list <- map(1:6, make_atac_plot)

ggsave("./sfigs/sfig9_atactracks.png", 
       cowplot::plot_grid(plotlist = atac_plot_list, ncol = 2),
       width = 6.5, height = 8)
