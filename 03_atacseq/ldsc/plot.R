library(tidyverse)
library(scico)

traits <- 
    read_tsv("./data/traits.txt", col_names = c("directory", "trait", "gwas", "ref")) |>
    select(gwas, trait)


results <- 
    read_tsv("./compiled_results.tsv") |>
    mutate(set = factor(set, levels = c("IL4", "TLR7", "BCR", "DN2")),
	   group = case_when(grepl("T2D|Height|HEIGHTz|Type_2_Diabetes|Schizophrenia|MDD|LDL|HDL|Covid19_Infection|Vaccination|cancer", gwas) ~ "control",
			     TRUE ~ "test")) |>
    left_join(traits, join_by(gwas)) |>
    arrange(trait, gwas, set) |>
    mutate(trait = fct_inorder(trait), gwas = fct_inorder(gwas))




p <- 
    ggplot(data = results, aes(x = set, y = gwas)) +
    geom_tile(aes(fill = tau_star)) +
    geom_text(data = filter(results, pfdr >= 0.01, pfdr <= 0.05),
	      aes(x = set, y = gwas, label = "*"), 
	      size = 10, fontface = "bold", size.unit = "pt", nudge_y = -.1) +
    geom_text(data = filter(results, pfdr <= 0.01),
	      aes(x = set, y = gwas, label = "**"), 
	      size = 10, fontface = "bold", size.unit = "pt", nudge_y = -.1) +
    scale_x_discrete(position = "top") +
    scale_fill_scico(palette = "vik", midpoint = 0) +
    facet_grid(rows = vars(group), scales = "free_y", space = "free_y") +
    theme_minimal() +
    theme(axis.title = element_blank(),
	  strip.text = element_blank(),
	  plot.background = element_rect(fill = "white", color = "white")
	  ) +
    labs(fill = "Tau*") +
    guides(fill = guide_colorbar(barheight = 12, barwidth = .5))

ggsave("./ldsc.png", p, width = 6)

