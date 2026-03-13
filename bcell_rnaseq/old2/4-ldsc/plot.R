library(tidyverse)

fig_colors <- 
    c("Day0" = "#BBBBBC", "TLR7" = "#52D464", "BCR" = "#4085D4", "DN2" = "#D44789")

ldsc_df <- 
    read_tsv("./compiled_results.tsv") |>
    mutate(gwas = fct_inorder(gwas),
	   set = factor(set, levels = c("Day0", "TLR7", "BCR", "DN2")))

ldsc_plot <- 
    ggplot(ldsc_df, aes(x = tau_hat, y = gwas)) +
    geom_col(aes(fill = set), position = "dodge") +
    geom_errorbar(aes(xmin = tau_hat - se*1.96, xmax = tau_hat + se*1.96),
		  width = .3, linewidth = .25) +
    geom_text(data = ldsc_df |> filter(pfdr <= 0.05),
	      aes(x = (tau_hat + se*1.96) * 1.1, y = gwas),
	      label = "*", nudge_y = -.2, fontface = "bold", color = "red") +
    scale_x_continuous(breaks = ~range(.),
		       labels = ~signif(., digits = 2)) +
    scale_fill_manual("Stim:", values = fig_colors) +
    facet_wrap(~set, nrow = 1) +
    theme_minimal() +
    theme(
	  axis.text.x = element_text(size = 8, hjust = c(0, 1)),
	  axis.ticks.x = element_line(linewidth = .25),
	  axis.ticks.length.x = unit(.1, "cm"),
	  axis.text.y = element_text(size = 8),
	  legend.position = "none",
	  plot.background = element_rect(fill = "white", color = "white")
	  ) +
    labs(y = NULL)

ggsave("./ldsc_ase_tauhat.png", ldsc_plot, height = 5)

ldsc_plot_norm <- 
    ggplot(ldsc_df, aes(x = tau_star, y = gwas)) +
    geom_col(aes(fill = set), position = "dodge") +
    scale_x_continuous(breaks = c(0, .1)) +
    scale_fill_manual("Stim:", values = fig_colors) +
    facet_wrap(~set, nrow = 1) +
    theme_minimal() +
    theme(
	  axis.text.x = element_text(size = 8),
	  axis.ticks.x = element_line(linewidth = .25),
	  axis.ticks.length.x = unit(.1, "cm"),
	  axis.text.y = element_text(size = 8),
	  legend.position = "none",
	  plot.background = element_rect(fill = "white", color = "white")
	  ) +
    labs(y = NULL)

ggsave("./ldsc_ase_taustar.png", ldsc_plot_norm, height = 5)

