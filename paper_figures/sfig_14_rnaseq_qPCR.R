library(DESeq2)
library(tidyverse)
library(ggbeeswarm)
library(patchwork)
library(ggh4x)

# Colors
stim_colors <- 
    "./figure_colors_final.txt" |>
    read_tsv() |>
    unite("stim", c(Condition, Time), sep = "_") |>
    select(stim, Hex) |>
    deframe()

stim_colors <- 
    stim_colors[c(1, 27:30)]

names(stim_colors)[1] <- "DN2c_0"

# Gene expression data
cpm_df <- 
    "../01_rnaseq_lowinput/2_dge/results/edger/cpm.tsv" |>
    read_tsv()


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

gene_set <- c("BLK", "IKZF2", "IRF5", "IRF7", "SNRPC", "WDFY4")

gene_merge <- 
    cpm_df |>
    filter(stim == "DN2", gene_name %in% gene_set) |>
    separate(sample_id, c("donor", "dummy", "timep"), sep = "_") |>
    select(donor, condition = stim, timep, gene = gene_name, value = obs_cpm) |>
    add_column(assay = "RNA-seq", .before = 1) |>
    mutate(
	   timep = str_remove(timep, "hrs"),
	   condition = recode(condition, "DN2" = "DN2c"),
	   condition = paste(condition, timep, sep = "_")
	   ) |>
    filter(timep %in% c("0", "24", "48", "72")) |>
    bind_rows(q_pcr_data) |>
    mutate(assay = fct_inorder(assay),
	   timep = factor(timep, levels = c("0", "24", "48", "72")))

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

ggsave("./sfigs/sfig14_rnaseq_qPCR.png", plot_merge, width = 6.5, height = 6.5)

#
#plot_list <- 
#    map(gene_set,
#        function(g) {
#            ggplot(data = gene_merge |> filter(assay == "RNA-seq", gene == g), 
#               aes(x = timep, y = value)) +
#            geom_quasirandom(aes(fill = condition),
#        		     method = "smiley", width = .2, 
#        		     shape = 21, stroke = .2, size = 2.5) +
#            scale_fill_manual(values = stim_colors) + 
#            scale_y_continuous(limits = c(0, NA),
#                               labels = ~str_pad(.x, width = 4, side = "left", pad = " ")
#            ) +
#            facet_wrap(~factor(gene, levels = gene_set),
#                       scale = "free_y", nrow = 1) +
#            theme_bw() +
#            theme(panel.grid.minor = element_blank(),
#        	  panel.grid.major.x = element_blank(),
#        	  axis.text = element_text(size = 8),
#        	  axis.title = element_text(size = 9),
#        	  strip.text = element_text(size = 9, angle = 0, face = "italic"),
#        	  legend.position = "none",
#        	  plot.background = element_rect(fill = "white", color = "white")) + 
#            labs(x = NULL, y = "Norm. counts")
#        }
#    )
#
#row1 <- plot_list[[1]] + plot_spacer() + plot_list[[2]] + plot_spacer() + plot_layout(nrow = 1, widths = c(1, 1.2, 1, 1.2))
#row2 <- plot_list[[3]] + plot_spacer() + plot_list[[4]] + plot_spacer() + plot_layout(nrow = 1, widths = c(1, 1.2, 1, 1.2))
#row3 <- plot_list[[5]] + plot_spacer() + plot_list[[6]] + plot_spacer() + plot_layout(nrow = 1, widths = c(1, 1.2, 1, 1.2))
#
#plot_out <- row1 / row2 / row3 + plot_layout(ncol = 1)
#
#ggsave("./supp_fig_RNAseq_selected.png", plot_out, width = 6.5, height = 6.5)
#
