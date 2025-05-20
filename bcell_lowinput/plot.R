library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(ggsci)
library(ggbeeswarm)
library(ComplexUpset)
library(patchwork)
library(gridExtra)

stim_colors <- 
    "../figure_colors.txt" |>
    read_tsv(col_names = c("condition", "color")) |>
    mutate(condition = paste0(condition, "hrs")) |>
    deframe()

all_stims <-
    names(stim_colors) |>
    str_split("_") |>
    map_chr(1) |>
    unique()

cpm_df <- 
    "./results/edger/cpm.tsv" |>
    read_tsv() |>
    group_by(sample_id, stim) |>
    nest() |>
    ungroup() |>
    separate(sample_id, c("sample_id", "dummy", "time"), sep = "_") |>
    mutate(hours = parse_number(time),
	   hours = factor(hours, levels = sort(unique(hours))),
           condition = paste(stim, time, sep = "_"),
	   stim = factor(stim, levels = all_stims)) |>
    unnest(cols = data) |>
    select(sample_id, condition, stim, hours, gene_id, gene_name, obs_cpm, obs_logcpm)

edger_results <- 
    read_tsv("./results/edger/results.tsv") |>
    mutate(stim = factor(stim, levels = all_stims))

selected_genes <- read_lines("./data/sle_curated_genes.txt")

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



selected_genes <- "IRF5"

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

stim_colors_0_names <- 
    all_stims[all_stims != "Unstim"] |>
    paste0("_0hrs")

stim_colors_0 <- rep(stim_colors[[1]], length(stim_colors_0_names))
names(stim_colors_0) <- stim_colors_0_names 


b_tfs_plot <-
    ggplot(data = cpm_plot_df |> filter(stim %in% c("CD40L", "DN2"))) +
    geom_quasirandom(aes(x = hours, y = obs_cpm, fill = condition),
		     method = "smiley", width = .15, 
		     shape = 21, stroke = .2, size = 2.75) +
#    geom_text(data = p_vals, 
#	      aes(x = 0.5, y = cpm * 1.25, label = p_lab),
#	      hjust = "inward", vjust = "inward", size = 4) +
    #scale_y_continuous(expand = expansion(mult = c(.2, .2))) +
    scale_fill_manual(values = c(stim_colors_0, stim_colors)) + 
    facet_grid(gene_name~fct_relevel(stim, levels(cpm_plot_df)),
	       scale = "free", space = "free_x") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
	  panel.grid.major.x = element_blank(),
	  panel.background = element_rect(fill = "grey96"),
	  strip.text.x = element_text(size = 10),
	  strip.text.y = element_text(angle = 0, size = 12, margin = margin(r = 0, l = 0)),
	  axis.text = element_text(size = 9),
	  axis.title = element_text(size = 12),
	  legend.position = "none") +
    labs(x = NULL, y = "Counts per million")

#ggsave("./plots/tbet_zeb2.png", b_tfs_plot, height = 6, width = 10)
ggsave("./plots/irf5.png", b_tfs_plot, height = 2, width = 4.5, dpi = 600)

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



# Differential expression
diff_expr <- 
    "./results/edger/diff_expr_all_times.tsv" |>
    read_tsv()

diff_expr_summ <-
    diff_expr |>
    group_by(group1, group2) |>
    summarise(n = sum(!is.na(gene_id))) |>
    ungroup() |>
    separate(group1, c("stim1", "t1"), sep = "\\.", remove = FALSE, convert = TRUE) |>
    separate(group2, c("stim2", "t2"), sep = "\\.", remove = FALSE, convert = TRUE) |>
    mutate_at(vars(stim1, stim2), ~recode(., "BCR_TLR7" = "BCR-TLR7")) |>
    mutate_at(vars(stim1, stim2), ~factor(., levels = all_stims)) |>
    complete(group1, group2, fill = list(n = NA)) |>
    arrange(stim1, stim2, t1, t2) |>
    mutate_at(vars(group1, group2), ~factor(., levels = unique(c(group1, group2))))

x_lines_df <- 
    diff_expr_summ |>
    filter(!is.na(n)) |>
    distinct(group1, stim1, t1) |>
    rowid_to_column() |>
    group_by(stim1) |>
    summarise(avg = mean(rowid), 
	      rowid = max(rowid) + 0.5) |>
    ungroup()

y_lines_df <- 
    diff_expr_summ |>
    filter(!is.na(n)) |>
    distinct(group2, stim2, t2) |>
    rowid_to_column() |>
    group_by(stim2) |>
    summarise(avg = mean(rowid), 
	      rowid = max(rowid) + 0.5) |>
    ungroup()

diff_plot <- 
    ggplot(data = diff_expr_summ, 
       aes(x = group1, y = group2)) +
    geom_tile(aes(fill = n), alpha = .9) +
    geom_vline(xintercept = head(x_lines_df$rowid, -1), 
	       color = "midnightblue", linewidth = .35) +
    geom_hline(yintercept = head(y_lines_df$rowid, -1), 
	       color = "midnightblue", linewidth = .35) +
    scale_x_discrete(labels = function(x) str_extract(x, "\\d+$")) +
    scale_y_discrete(labels = function(x) str_extract(x, "\\d+$")) +
    scale_fill_viridis_c(option = "cividis",
			 na.value = "white",
			 labels = scales::comma) +
    theme_minimal() +
    theme(
	  axis.title = element_blank(),
	  panel.grid.major = element_line(color = "black", linewidth = .35),
	  legend.position = "top",
	  legend.title.position = "top",
	  plot.margin = margin(b = 0.25, l = 0.75, unit = "in"),
	  plot.background = element_rect(fill = "white", color = "white")) +
    labs(fill = "Number of DE genes:") +
    guides(fill = guide_colorbar(barwidth = 15, barheight = .5)) +
    annotate("text", x = x_lines_df$avg, y = -1, label = x_lines_df$stim1, 
	     size = 9, size.unit = "pt") +
    annotate("text", x = -1, y = y_lines_df$avg, label = y_lines_df$stim2, 
	     size = 9, size.unit = "pt", hjust = 1) +
    coord_cartesian(xlim = c(1, 28), ylim = c(1, 28), clip = "off")

ggsave("./paper_plots/", diff_plot, width = 6, height = 6) 


## splines
diff_expr_splines <- 
    "./results/edger/diff_expr_splines.tsv" |>
    read_tsv()

diff_expr_splines_summ <- 
    diff_expr_splines |>
    filter(!is.na(PValue)) |>    
    count(group1, group2) |>
    bind_rows(filter(diff_expr_splines, is.na(PValue)) |> mutate(n = 0)) |>
    mutate_at(vars(group1, group2), ~recode(., "BCR_TLR7" = "BCR-TLR7")) |>
    mutate_at(vars(group1, group2), ~factor(., levels = all_stims)) |>
    arrange(group1, group2)

testp <- 
    ggplot(data = diff_expr_splines_summ,
	   aes(x = group1, y = group2)) +
    geom_tile(aes(fill = n)) +
    geom_text(aes(label = n), color = "white", fontface = "bold") +
    scale_fill_viridis_c(option = "cividis",
			 na.value = "white",
			 labels = scales::comma) +
    theme_minimal() +
    theme(axis.title = element_blank(),
	  plot.background = element_rect(fill = "white", color = "white")) +
    guides(fill = "none")
    
ggsave("./plots/testp.png", testp, width = 4, height = 4)

# RPS24
salmon_tx <-
    read_tsv("./data/expression_pooled_reps_transcriptlevel.tsv") |>
    
rps24 <- 
    salmon_tx |>
    filter(gene_name == "RPS24") |>
    separate(sample_id, c("donor_id", "stim", "timep"), sep = "_") |>
    mutate(timep = parse_number(timep),
	   stim = factor(stim, levels = c("Unstim", "IL4", "CD40L", "TLR7", "TLR9", "BCR", "BCR-TLR7", "DN2")))

rps24_gene <- 
    rps24 |>
    group_by(donor_id, stim, timep, gene_id, gene_name) |>
    summarise(tpm = sum(tpm)) |>
    ungroup()

rps_plot <- 
    ggplot(rps24_gene |> filter(stim == "DN2"), 
	   aes(x = factor(timep), y = tpm)) +
    geom_quasirandom(method = "smiley", width = .2) +
    facet_wrap(~gene_name, ncol = 2) +
    theme_bw() +
    theme(panel.grid.minor.y = element_blank())


rps_plot_tx <- 
    rps24 |> 
    filter(stim == "DN2") |> 
    group_by(tx_id) |>
    filter(sum(tpm >= 5) > 2) |>
    ungroup() |>
    ggplot(aes(x = factor(timep), y = tpm)) +
    geom_quasirandom(method = "smiley", width = .2) +
    facet_wrap(~tx_id, ncol = 2, scales = "free_y") +
    theme_bw() +
    theme(panel.grid.minor.y = element_blank())

rps_out <- 
    plot_grid(plot_grid(rps_plot, NULL, ncol = 1, rel_heights = c(.64, 1)),
	  rps_plot_tx, 
	  nrow = 1, rel_widths = c(.5, 1)) +
    theme(plot.background = element_rect(fill = "white", color = "white"))

ggsave("./plots/rps24_dn2.png", rps_out, width = 7.5, height = 5)

