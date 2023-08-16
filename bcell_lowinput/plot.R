library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(ggsci)
library(ggbeeswarm)
library(patchwork)

stim_colors <- 
    c("CD40L_0hrs" = "#dbdbdb", 
     "TLR9_0hrs" = "#dbdbdb", 
     "TLR7_0hrs" = "#dbdbdb",
     "BCR_0hrs" = "#dbdbdb",
     "BCR+TLR7_0hrs" = "#dbdbdb",
     "DN2_0hrs" = "#dbdbdb",
     "CD40L_4hrs" = "#ffeba6",
     "CD40L_24hrs" = "#ffde68",
     "CD40L_48hrs" = "goldenrod3",
     "CD40L_72hrs" = "goldenrod4",
     "BCR_4hrs" = "#a1c2ed",
     "BCR_24hrs" = "#6996e3",
     "BCR_48hrs" = "#4060c8",
     "BCR_72hrs" = "#0404bf",
     "TLR7_4hrs" = "#98ab76",
     "TLR7_24hrs" = "#748f46",
     "TLR7_48hrs" = "#47632a",
     "TLR7_72hrs" = "#275024",
     "TLR9_4hrs" = "#a876d9",
     "TLR9_24hrs" = "#955bd0",
     "TLR9_48hrs" = "#803ec8",
     "TLR9_72hrs" = "#691dbf",
     "BCR+TLR7_4hrs" = "#ffa7db",
     "BCR+TLR7_24hrs" = "#ff86d0",
     "BCR+TLR7_48hrs" = "#ff61c4",
     "BCR+TLR7_72hrs" = "#ff2bb8",
     "DN2_4hrs" = "#e6907a",
     "DN2_24hrs" = "#d76b51",
     "DN2_48hrs" = "#c5432a",
     "DN2_72hrs" = "#b00000")

cpm_df <- 
    "./results/edger/cpm.tsv" |>
    read_tsv() |>
    group_by(sample_id, stim) |>
    nest() |>
    ungroup() |>
    separate(sample_id, c("sample_id", "dummy", "time"), sep = "_") |>
    mutate(hours = parse_number(time),
	   hours = factor(hours, levels = sort(unique(hours))),
	   stim = recode(stim, "BCR-TLR7" = "BCR+TLR7"),
           condition = paste(stim, time, sep = "_"),
	   stim = factor(stim, levels = c("CD40L", "TLR7", "TLR9", "BCR", "BCR+TLR7", "DN2"))) |>
    unnest(cols = data) |>
    select(sample_id, condition, stim, hours, gene_id, gene_name, obs_cpm, fit_cpm)

sle_genes <- read_tsv("../bcell_scrna/reported_genes.tsv")

edger_results <- 
    read_tsv("./results/edger/results.tsv") |>
    mutate(stim = recode(stim, "BCR-TLR7" = "BCR+TLR7"),
	   stim = factor(stim, levels = c("CD40L", "TLR7", "TLR9", "BCR", "BCR+TLR7", "DN2")))

selected_genes <- 
    c("BACH2", "BANK1", "IRF5", "IRF8", "JAZF1", "PTPN22", "PXK", 
      "SOCS1", "STAT1", "TNFSF4", "WDFY4", "IL10")

cpm_plot_df <- cpm_df |>
    filter(gene_name %in% selected_genes)

p_vals <- edger_results |>
    inner_join(distinct(cpm_plot_df, stim, gene_id, gene_name)) |>
    select(gene_id, gene_name, stim, p = PValue, fdr = FDR) |>
    arrange(gene_name, stim) |>
    left_join(summarise(cpm_plot_df, cpm = max(obs_cpm), .by = c(gene_id, gene_name))) |>
    mutate(p_lab = format(p, format = "e", digits = 2),
	   p_lab = paste("p =", p_lab),
	   p_lab = ifelse(fdr < 0.05, paste(p_lab, "*"), p_lab))

out_plot <-
    cpm_plot_df |>
    group_by(gene_id, gene_name) |>
    nest() |>
    ungroup() |>
    mutate(s = ntile(gene_name, 2)) |>
    unnest(cols = data) |>
    group_split(s) |>
    map(function(x) ggplot(data = x) +
	geom_quasirandom(aes(x = hours, y = obs_cpm, fill = condition),
			 method = "smiley", width = .2, 
			 shape = 21, stroke = .2, size = 4) +
	geom_text(data = p_vals |> filter(gene_name %in% x$gene_name), 
		  aes(x = 0.5, y = cpm * 1.25, label = p_lab),
		  hjust = "inward", vjust = "inward", size = 4) +
	scale_fill_manual(values = stim_colors) + 
	facet_grid(gene_name~fct_relevel(stim, levels(cpm_plot_df)),
		   scale = "free_y") +
	theme(panel.grid.minor = element_blank(),
	      panel.grid.major.x = element_blank(),
	      panel.background = element_rect(fill = "grey96"),
	      strip.text.x = element_text(size = 12),
	      strip.text.y = element_text(angle = 0, size = 12),
	      axis.text = element_text(size = 12),
	      axis.title = element_text(size = 14),
	      legend.position = "none") +
	labs(x = NULL, y = "Normalized counts"))

ggsave("./plots/timecourse_topgenes.png",
       out_plot[[1]] + out_plot[[2]] + plot_layout(widths = c(1, 1)),
       width = 19, height = 7, dpi = 600)



# T-bet and ZEB2
selected_genes <- c("TBX21", "ZEB2") 

cpm_plot_df <- cpm_df |>
    filter(gene_name %in% selected_genes)

p_vals <- edger_results |>
    inner_join(distinct(cpm_plot_df, stim, gene_id, gene_name)) |>
    select(gene_id, gene_name, stim, p = PValue, fdr = FDR) |>
    arrange(gene_name, stim) |>
    left_join(summarise(cpm_plot_df, cpm = max(obs_cpm), .by = c(gene_id, gene_name))) |>
    mutate(p_lab = format(p, format = "e", digits = 2),
	   p_lab = paste("p =", p_lab),
	   p_lab = ifelse(fdr < 0.05, paste(p_lab, "*"), p_lab))

b_tfs_plot <-
    ggplot(data = cpm_plot_df) +
    geom_quasirandom(aes(x = hours, y = obs_cpm, fill = condition),
		     method = "smiley", width = .2, 
		     shape = 21, stroke = .2, size = 3.5) +
    geom_text(data = p_vals, 
	      aes(x = 0.5, y = cpm * 1.25, label = p_lab),
	      hjust = "inward", vjust = "inward", size = 4) +
    scale_fill_manual(values = stim_colors) + 
    facet_grid(gene_name~fct_relevel(stim, levels(cpm_plot_df)),
	       scale = "free_y") +
    theme(panel.grid.minor = element_blank(),
	  panel.grid.major.x = element_blank(),
	  panel.background = element_rect(fill = "grey96"),
	  strip.text.x = element_text(size = 12),
	  strip.text.y = element_text(angle = 0, size = 12),
	  axis.text = element_text(size = 12),
	  axis.title = element_text(size = 14),
	  legend.position = "none") +
    labs(x = NULL, y = "Normalized counts")

ggsave("./plots/tbet_zeb2.png", b_tfs_plot, height = 4, width = 10)


