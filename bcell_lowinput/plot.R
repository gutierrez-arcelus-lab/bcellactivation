library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(ggsci)
library(ggbeeswarm)
library(ComplexUpset)
library(patchwork)
library(gridExtra)

stim_colors <- 
    c("IL4_0hrs" = "#dbdbdb",
      "CD40L_0hrs" = "#dbdbdb", 
      "TLR9_0hrs" = "#dbdbdb", 
      "TLR7_0hrs" = "#dbdbdb",
      "BCR_0hrs" = "#dbdbdb",
      "BCR+TLR7_0hrs" = "#dbdbdb",
      "DN2_0hrs" = "#dbdbdb",
      "IL4_4hrs" = "grey50",
      "IL4_24hrs" = "black",
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

all_stims <- c("IL4", "CD40L", "TLR7", "TLR9", "BCR", "BCR+TLR7", "DN2")

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
	   stim = factor(stim, levels = all_stims)) |>
    unnest(cols = data) |>
    select(sample_id, condition, stim, hours, gene_id, gene_name, obs_cpm, obs_logcpm)

edger_results <- 
    read_tsv("./results/edger/results.tsv") |>
    mutate(stim = recode(stim, "BCR-TLR7" = "BCR+TLR7"),
	   stim = factor(stim, levels = all_stims))

selected_genes <-
    c(
      "ARID5B",
      "ATG5",
      "BANK1",
      "BLK", 
      "CD44",
      "CSK",
      "DHCR7",
      "ETS1",
      "FCGR2A",
      "LYST",
      "IFIH1",
      "IKBKE",
      "IKZF1",
      "IKZF2",
      "IKZF3",
      "IL10",
      "IL12A",
      "IRAK1",
      "IRF5",
      "IRF7",
      "IRF8",
      "ITGAM",
      "ITGAX",
      "JAZF1",
      "MECP2",
      "MIR146A",
      "MIR3142HG",
      "NCF2",
      "PTPN22",
      "PRDM1",
      "PXK",
      "RAD51B",
      "SH2B3",
      "SLC15A4",
      "SOCS1",
      "SPRED2",
      "STAT1", 
      "STAT4",
      "TASL",
      "TCF7",
      "TNFAIP3",
      "TNFSF4",
      "TNIP1",
      "TREX1",
      "TYK2",
      "UHRF1BP1",
      "UBE2L3",
      "WDFY4"
    )

cpm_plot_df <- cpm_df |>
    filter(gene_name %in% selected_genes)

p_vals <- edger_results |>
    inner_join(distinct(cpm_plot_df, stim, gene_id, gene_name)) |>
    select(gene_id, gene_name, stim, p = PValue, fdr = FDR) |>
    arrange(gene_name, stim) |>
    left_join(summarise(cpm_plot_df, cpm = max(obs_cpm), .by = c(gene_id, gene_name))) |>
    mutate(p_lab = format(p, format = "e", digits = 2),
	   p_lab = ifelse(fdr < 0.05, paste(p_lab, "*"), p_lab))

out_plot <-
    cpm_plot_df |>
    {function(x) split(x, x$gene_name)}() |>
    map(~ggplot(data = .x) +
	geom_quasirandom(aes(x = hours, y = obs_cpm, fill = condition),
			 method = "smiley", width = .2, 
			 shape = 21, stroke = .2, size = 3.5) +
	geom_text(data = filter(p_vals, gene_name %in% unique(.x$gene_name)),
		  aes(x = 0.5, y = cpm * 1.25, label = p_lab),
		  hjust = "inward", vjust = "inward", size = 4) +
	scale_fill_manual(values = stim_colors) + 
	facet_grid(.~fct_relevel(stim, levels(cpm_plot_df)),
		   scale = "free", space = "free_x", drop = TRUE) +
	theme(panel.grid.minor = element_blank(),
	      panel.grid.major.x = element_blank(),
	      panel.background = element_rect(fill = "grey96"),
	      strip.text.x = element_text(size = 12),
	      strip.text.y = element_text(angle = 0, size = 12),
	      axis.text = element_text(size = 12),
	      axis.title = element_text(size = 14),
	      legend.position = "none",
	      plot.margin = margin(3.2, 0.1, 3.2, 0.1, unit = "in")) +
	labs(x = NULL, y = "Normalized counts", title = unique(.x$gene_name)))

walk(names(out_plot), 
     ~ggsave(sprintf("./plots/plots_timecourse/%s.pdf", .x), 
	     out_plot[[.x]], 
	     width = 11,
	     height = 8.5,
	     dpi = 300))



# T-bet and ZEB2
selected_genes <- c("TBX21", "ZEB2", "SLC15A4") 

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

ggsave("./plots/tbet_zeb2.png", b_tfs_plot, height = 6, width = 10)

# Test genes
cpm_df <- 
    "./results/edger/cpm_2.tsv" |>
    read_tsv() |>
    group_by(sample_id, stim) |>
    nest() |>
    ungroup() |>
    separate(sample_id, c("sample_id", "dummy", "time"), sep = "_") |>
    mutate(hours = parse_number(time),
	   hours = factor(hours, levels = sort(unique(hours))),
	   stim = recode(stim, "BCR-TLR7" = "BCR+TLR7"),
           condition = paste(stim, time, sep = "_"),
	   stim = factor(stim, levels = all_stims)) |>
    unnest(cols = data) |>
    select(sample_id, condition, stim, hours, gene_id, gene_name, obs_logcpm)

edger_results <- 
    read_tsv("./results/edger/results_2.tsv") |>
    mutate(stim = recode(stim, "BCR-TLR7" = "BCR+TLR7"),
	   stim = factor(stim, levels = all_stims))

selected_genes <- c("MAPK10") 

cpm_plot_df <- cpm_df |>
    filter(gene_name %in% selected_genes)

p_vals <- edger_results |>
    inner_join(distinct(cpm_plot_df, stim, gene_id, gene_name)) |>
    select(gene_id, gene_name, stim, p = PValue, fdr = FDR) |>
    arrange(gene_name, stim) |>
    left_join(summarise(cpm_plot_df, cpm = max(obs_logcpm), .by = c(gene_id, gene_name))) |>
    mutate(p_lab = format(p, format = "e", digits = 2),
	   p_lab = paste("p =", p_lab),
	   p_lab = ifelse(fdr < 0.05, paste(p_lab, "*"), p_lab))

b_tfs_plot <-
    ggplot(data = cpm_plot_df) +
    geom_quasirandom(aes(x = hours, y = obs_logcpm, fill = condition),
		     method = "smiley", width = .2, 
		     shape = 21, stroke = .2, size = 3.5) +
    geom_text(data = p_vals, 
	      aes(x = 0.5, y = cpm * 1.25, label = p_lab),
	      hjust = "inward", vjust = "inward", size = 4) +
    scale_fill_manual(values = stim_colors) + 
    facet_grid(gene_name~fct_relevel(stim, levels(cpm_plot_df)),
	       scale = "free") +
    theme(panel.grid.minor = element_blank(),
	  panel.grid.major.x = element_blank(),
	  panel.background = element_rect(fill = "grey96"),
	  strip.text.x = element_text(size = 12),
	  strip.text.y = element_text(angle = 0, size = 12),
	  axis.text = element_text(size = 12),
	  axis.title = element_text(size = 14),
	  legend.position = "none") +
    labs(x = NULL, y = "Normalized counts")

ggsave("./plots/test_2.png", b_tfs_plot, height = 2.5, width = 12)




# UpSet plot
edger_memb <- edger_results |>
    filter(stim != "IL4") |>
    select(stim, gene_id, FDR) |>
    pivot_wider(names_from = stim, values_from = FDR) |>
    mutate_at(vars(-gene_id), ~ifelse(. < 0.05, 1, 0))

p <- 
    upset(edger_memb, names(edger_memb)[-1], 
	  base_annotations = list('Intersection size' = intersection_size(text = list(vjust = 0,
										      hjust = 0.5,
										      angle = 0,
										      size = 2.5)) +
				  theme(panel.grid = element_blank())),
	  set_sizes = 
	      (upset_set_size() + 
	       theme(axis.text.x = element_text(angle = 90),
		     axis.ticks.x = element_line()) +
	       geom_text(aes(label = after_stat(count)), 
			 stat = "count", 
			 size = 3, hjust = 0, color = "white")),
	  name = NULL, width_ratio = 0.125, min_size = 100)

ggsave("./plots/upset.png", p, width = 10, height = 5)
