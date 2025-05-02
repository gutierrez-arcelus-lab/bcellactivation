library(tidyverse)
library(glue)
library(cowplot)

stim_i <- commandArgs(TRUE)[1]

module_sizes <-
    glue("../bcell_lowinput/wgcna/data/{stim_i}_modules.tsv") |>
    read_tsv() |>
    count(module) |>
    arrange(desc(n)) |>
    filter(module != "grey") |>
    rowid_to_column("ix") |>
    select(module, ix) |>
    mutate(ix = paste("Module", ix)) |>
    mutate(ix = fct_inorder(ix))

eigengenes_df <-
    glue("../bcell_lowinput/wgcna/data/{stim_i}_eigen.tsv") |>
    read_tsv() |>
    select(-MEgrey) |>
    separate(sample_name, 
	     c("donor_id", "stim", "time"), 
	     sep = "_") |>
    pivot_longer(starts_with("ME"), names_to = "module") |>
    mutate(module = str_remove(module, "^ME"),
	   time = parse_number(time),
	   time = factor(time, levels = sort(unique(time)))) |>
    left_join(module_sizes, join_by(module)) |>
    select(donor_id, module_ix = ix, module, time, value) |>
    arrange(module_ix, donor_id, time)

modules_df <-
    eigengenes_df |>
    group_by(module, module_ix, time) |>
    summarise(value = mean(value)) |>
    ungroup() |>
    arrange(module_ix, time)

kim_df <- 
    glue("../bcell_lowinput/wgcna/data/{stim_i}_kim.tsv") |>
    read_tsv() |>
    left_join(module_sizes, join_by(module)) |> 
    select(gene_id, gene_name, module, module_ix = ix, kim) |>
    arrange(module_ix, desc(kim)) |>
    group_by(module) |>
    top_n(5, kim) |>
    ungroup()

module_plot <-
    ggplot(modules_df) +
    geom_line(aes(x = time, y = value, group = module),
	      linewidth = 1.25) +
    facet_wrap(~module_ix, ncol = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 8),
	  axis.text.y = element_text(size = 8),
	  axis.title.x = element_blank(),
	  axis.title.y = element_text(size = 9),
	  panel.grid.major.x = element_line(linetype = 2, color = "black", linewidth = .1),
	  panel.grid.major.y = element_blank(),
	  panel.grid.minor.y = element_blank(),
	  panel.spacing.y = unit(0.1, "lines"),
	  strip.text = element_text(size = 8, margin = margin(t = 0, b = 0)),
	  strip.clip = "off",
	  plot.margin = margin(0, 0, 0, 0, unit = "lines")) +
    labs(y = "Eigengene expression")

kim_facet_labels <- 
    rep(" ", n_distinct(kim_df$module_ix)) |>
    setNames(as.character(unique(kim_df$module_ix)))

kim_plot <-
    ggplot(data = kim_df |> 
	   arrange(module_ix, kim) |>
	   mutate(gene_name = fct_inorder(gene_name)),
	   aes(x = " ", y = gene_name)) +
    geom_point(aes(fill = kim), shape = 21, size = 2.5, show.legend = FALSE) +
    scale_x_discrete(drop = FALSE, expand = c(0, 0), limits = c(" ", " ")) +
    scale_y_discrete(position = "right", expand = c(0, 0)) +
    scale_fill_continuous(low = "beige", high = "firebrick") +
    facet_wrap(~module_ix, scale = "free_y", ncol = 1,
	       labeller = as_labeller(kim_facet_labels)) + 
    theme_minimal() +
    theme(axis.text.y.right = element_text(size = 7, face = "italic", 
					   margin = margin(l = .5, unit = "lines")),
	  panel.grid = element_blank(),
	  panel.spacing.y = unit(0.1, "lines"),
	  plot.margin = margin(0, 0, 0, 0, unit = "lines"),
	  strip.text = element_text(size = 8.5, margin = margin(t = 0, b = 0)),
	  strip.clip = "off"
	  ) +
    labs(x = NULL, y = NULL) +
    coord_cartesian(clip = "off")

module_grid <- plot_grid(module_plot, NULL, kim_plot, nrow = 1, rel_widths = c(1, .1, .5))

go_res <-
    glue("../bcell_lowinput/wgcna/data/{stim_i}_go.tsv") |>
    read_tsv() |>
    group_by(module) |>
    top_n(5, -log10(pvalue)) |>
    ungroup() |>
    mutate(Description = str_trunc(Description, width = 36),
	   gene_r = map_dbl(GeneRatio, ~eval(parse(text = .)))) |>
    right_join(module_sizes, join_by(module)) |>
    select(module = ix, Description, pvalue, gene_r) |>
    arrange(module, pvalue) |>
    mutate(Description = fct_inorder(Description))

go_plot <- 
    ggplot(data = go_res, 
       aes(x = " ", y = Description)) +
    geom_point(aes(fill = -log10(pvalue), size = gene_r), 
	       shape = 21, stroke = .5) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(position = "right") +
    scale_fill_gradient(low = "white", high = "black") +
    scale_size(range = c(.5, 2.5)) +
    ggforce::facet_col(vars(module), scales = "free_y", drop = FALSE) + 
    theme_minimal() +
    theme(
	  axis.text.y.right = element_text(size = 7, margin = margin(r = .5, l = 0, unit = "lines")),
	  axis.title = element_blank(),
	  legend.text = element_text(size = 8),
	  legend.title = element_text(size = 8),
	  legend.key.size = unit(.25, "cm"),
	  legend.margin = margin(r = 0.1, l = 0.1, unit = "lines"),
	  panel.grid = element_blank(),
	  plot.margin = margin(0, 0, 0, 0, unit = "lines"),
	  strip.clip = "off",
	  strip.text = element_text(size = 8, margin = margin(t = 0, b = 0, l = 2, unit = "lines")),
	  panel.spacing.y = unit(0.1, "lines")
	  ) +
    guides(fill = guide_colorbar("logP:", barheight = 5, barwidth = .25),
	   size = guide_legend("Gene\nRatio:"))


ggsave(glue("./sfig_modules_{stim_i}.png"),
       plot_grid(module_grid, go_plot, nrow = 1, rel_widths = c(.8, 1)) + 
	   theme(plot.background = element_rect(fill = "white", color = "white")),
       width = 5.5, height = 5)
	    
