library(tidyverse)
#library(glue)
library(readxl)
#library(cowplot)
#library(edgeR)
#library(fgsea)
#library(mirai)

dis_state <- 
    read_excel("./data/mmc2.xlsx", sheet = 2, skip = 2) |>
    janitor::clean_names()

dis_activ <- 
    read_excel("./data/mmc2.xlsx", sheet = 3, skip = 2) |>
    janitor::clean_names()

tf_enrich <- 
    c("BATF", "ATF3", "FOSL1", "JUN", "JUNB", "FOSL2", "BACH2", "MAFK", "BACH1", "NFE2L2",
      "STAT6", "STAT5A", "STAT5B", "STAT1", "STAT4", "STAT3", 
      "REL", "RELA", "NFKB1", "NFKB2", "NFATC1", "NFATC2", 
      "PRDM1", "EGR1", "EGR2", "BCL6", "ZBTB33",
      "RUNX1", 
      "EBF1",
      "TBX21",
      "NR4A1",
      "SPI1",
      "POU2F2", "POU5F1",
      "ZNF143", "RBPJ",
      "MEF2A", "MEF2C")

cell_types_b <- c("Naive B", "USM B", "SM B", "DN B", "Plasmablast")

tf_nakano_data <- 
    bind_rows("Disease state" = dis_state, "Disease activity" = dis_activ, .id = "contrast") |> 
    filter(cell_type %in% cell_types_b, gene %in% tf_enrich) |>
    mutate(cell_type = factor(cell_type, levels = cell_types_b),
	   gene = factor(gene, levels = rev(tf_enrich)))

tf_state_plot <- 
    ggplot(tf_nakano_data) +
    geom_tile(aes(x = cell_type, y = gene, fill = log_fc)) +
    geom_text(data = filter(tf_nakano_data, fdr <= 0.05),
	      aes(x = cell_type, y = gene, label = "*"), 
	      vjust = .75) +
    scale_fill_gradient2(
			 low = "#2166AC",    # blue
			 mid = "white", 
			 high = "#B2182B",   # red
			 midpoint = 0,
			 limits = c(-2.5, 2.5),
			 oob = scales::squish
    ) +
    facet_wrap(~fct_inorder(contrast), nrow = 1) +
    theme_minimal() +
    theme(axis.text = element_text(size = 7),
	  axis.title = element_text(size = 8),
	  plot.background = element_rect(fill = "white", color = "white")) +
    labs(x = "Cell type in Nakano et al.", y = NULL) +
    guides(fill = guide_colorbar(barwidth = .5))

ggsave("./sfigs/sfig11_tf_nakano_sigs.png", tf_state_plot, width = 6.5, height = 4)


