library(DESeq2)
library(tidyverse)
library(ggbeeswarm)
library(patchwork)
library(ggh4x)

# Colors
stim_colors <- 
    "./figure_colors.txt" |>
    read_tsv(col_names = c("stim", "time", "color")) |>
    unite("condition", c(stim, time), sep = "_") |>
    deframe()

stim_colors_0_names <- 
    keep(names(stim_colors), !grepl("Unstim", names(stim_colors))) |>
    str_remove("_\\d+hrs$") |>
    unique()

stim_colors_0 <- rep(stim_colors[[1]], length(stim_colors_0_names))
names(stim_colors_0) <- stim_colors_0_names

stim_colors <- c(stim_colors, stim_colors_0)

# Gene expression data
load("../../shinybcells/data/gene_exp.rda")

# PCR
pcr_colors <- 
    c(stim_colors[1], 
      "DN2c_24_pcr" = "#9dc6e0",
      "DN2c_72_pcr" = "#004c6d", 
      "DN2c_48_pcr" = "#5886a5")

q_pcr_data <- 
    readxl::read_excel("./data/Gene kinetics.xlsx") |>
    janitor::clean_names() |>
    select(gene = x1, everything()) |>
    pivot_longer(-gene) |>
    mutate(name = case_when(grepl("^x\\d+$", name) ~ NA,
                            TRUE ~ name)) |>
    fill(name) |>
    mutate(timep = fct_inorder(as.character(parse_number(name)))) |>
    group_by(gene, name) |>
    mutate(donor = paste0("PCR", seq_len(n()))) |>
    ungroup() |>
    add_column(assay = "qPCR", .before = 1) |>
    mutate(condition = case_when(timep == "0" ~ "Unstim_0_pcr",
                                 timep == "24" ~ "DN2c_24_pcr",
                                 timep == "48" ~ "DN2c_48_pcr",
                                 timep == "72" ~ "DN2c_72_pcr")) |>
    select(assay, donor, condition, timep, gene, value)
                                 

#
gene_set <- c("BLK", "IKZF2", "IRF5", "IRF7", "SNRPC", "WDFY4")

gene_exp_select <- 
    gene_exp |>
    select(donor, condition, stim, timep, all_of(gene_set)) |>
    pivot_longer(-(donor:timep), names_to = "gene")

norm_factors <- 
    gene_exp_select |>
    filter(condition == "Unstim_0") |>
    select(donor, gene, base_value = value)

gene_exp_norm <-
    gene_exp_select |>
    filter(condition == "Unstim_0" | stim == "DN2c") |>
    left_join(norm_factors, join_by(donor, gene)) |>
    mutate(norm_exp = value/base_value)

gene_merge <- 
    select(gene_exp_norm, donor, condition, timep, gene, value) |>
    add_column(assay = "RNA-seq", .before = 1) |>
    filter(timep %in% c("0", "24", "48", "72")) |>
    bind_rows(q_pcr_data) |>
    mutate(assay = fct_inorder(assay))

plot_merge <- 
    ggplot(data = gene_merge, 
       aes(x = timep, y = value)) +
    geom_quasirandom(aes(fill = condition),
		     method = "smiley", width = .2, 
		     shape = 21, stroke = .2, size = 2.5) +
    scale_fill_manual(values = c(stim_colors, pcr_colors)) + 
    scale_y_continuous(limits = c(0, NA),
                       labels = ~str_pad(.x, width = 4, side = "left", pad = " ")
    ) +
    facet_nested_wrap(vars(gene, assay), ncol = 4, scales = "free") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
	  panel.grid.major.x = element_blank(),
	  axis.text = element_text(size = 8),
	  axis.title = element_text(size = 9),
	  strip.text = element_text(size = 9, angle = 0, face = "italic"),
	  legend.position = "none",
	  plot.background = element_rect(fill = "white", color = "white")) + 
    labs(x = NULL, y = "Expression")

ggsave("./supp_fig_RNAseq_selected.png", plot_merge, width = 6.5, height = 6.5)


plot_list <- 
    map(gene_set,
        function(g) {
            ggplot(data = gene_merge |> filter(assay == "RNA-seq", gene == g), 
               aes(x = timep, y = value)) +
            geom_quasirandom(aes(fill = condition),
        		     method = "smiley", width = .2, 
        		     shape = 21, stroke = .2, size = 2.5) +
            scale_fill_manual(values = stim_colors) + 
            scale_y_continuous(limits = c(0, NA),
                               labels = ~str_pad(.x, width = 4, side = "left", pad = " ")
            ) +
            facet_wrap(~factor(gene, levels = gene_set),
                       scale = "free_y", nrow = 1) +
            theme_bw() +
            theme(panel.grid.minor = element_blank(),
        	  panel.grid.major.x = element_blank(),
        	  axis.text = element_text(size = 8),
        	  axis.title = element_text(size = 9),
        	  strip.text = element_text(size = 9, angle = 0, face = "italic"),
        	  legend.position = "none",
        	  plot.background = element_rect(fill = "white", color = "white")) + 
            labs(x = NULL, y = "Norm. counts")
        }
    )

row1 <- plot_list[[1]] + plot_spacer() + plot_list[[2]] + plot_spacer() + plot_layout(nrow = 1, widths = c(1, 1.2, 1, 1.2))
row2 <- plot_list[[3]] + plot_spacer() + plot_list[[4]] + plot_spacer() + plot_layout(nrow = 1, widths = c(1, 1.2, 1, 1.2))
row3 <- plot_list[[5]] + plot_spacer() + plot_list[[6]] + plot_spacer() + plot_layout(nrow = 1, widths = c(1, 1.2, 1, 1.2))

plot_out <- row1 / row2 / row3 + plot_layout(ncol = 1)

ggsave("./supp_fig_RNAseq_selected.png", plot_out, width = 6.5, height = 6.5)

# ATAC-seq
library(AnnotationHub)
library(locuszoomr)

select <- dplyr::select
filter <- dplyr::filter

ah <- AnnotationHub()
ens_data <- ah[["AH98047"]]

bigwigs <- 
    list.files("../../shinybcells/inst/extdata/",
	       pattern = "*_filtered.bigWig",
	       full.names = TRUE)
        
names(bigwigs) <- sub("^([^_]+_\\d+).+$", "\\1", basename(bigwigs))

load("../../shinybcells/data/atac_colors.rda")

da_files <- 
    list.files("../atacseq/results_deseq2/mRp",
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

ggsave("./test_atac_plot.png", 
       cowplot::plot_grid(plotlist = atac_plot_list, ncol = 2),
       width = 6.5, height = 8)
