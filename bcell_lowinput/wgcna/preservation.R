library(tidyverse)
library(WGCNA)

multi_expr <-
    list(
	 "CD40L"    = list(data = read_rds("./data/CD40L_counts.rds")),
	 "TLR9"     = list(data = read_rds("./data/TLR9_counts.rds")),
	 "TLR7"     = list(data = read_rds("./data/TLR7_counts.rds")),
	 "BCR"      = list(data = read_rds("./data/BCR_counts.rds")),
	 "BCR-TLR7" = list(data = read_rds("./data/BCR-TLR7_counts.rds")),
	 "DN2"      = list(data = read_rds("./data/DN2_counts.rds"))
    )

multi_color <- 
    list(
	 "CD40L"    = read_tsv("./data/CD40L_modules.tsv") |> deframe(),
	 "TLR9"     = read_tsv("./data/TLR9_modules.tsv") |> deframe(),
	 "TLR7"     = read_tsv("./data/TLR7_modules.tsv") |> deframe(),
	 "BCR"      = read_tsv("./data/BCR_modules.tsv") |> deframe(),
	 "BCR-TLR7" = read_tsv("./data/BCR-TLR7_modules.tsv") |> deframe(),
	 "DN2"      = read_tsv("./data/DN2_modules.tsv") |> deframe()
    )

mp <- 
    modulePreservation(multi_expr, multi_color, 
		       referenceNetworks = 1:6,
		       loadPermutedStatistics = FALSE,
		       networkType = "signed",
		       nPermutations = 100,
		       randomSeed = 1,
		       verbose = 3,
		       parallelCalculation = TRUE
    )

write_rds(mp, "data/module_preservation.rds")

mp <- read_rds("data/module_preservation.rds")

z_tidy <- 
    mp$preservation$Z |>
    list_flatten() |>
    map_dfr(~as_tibble(., rownames = "module"), .id = "pair") |>
    select(pair, module, ngenes = moduleSize, zsumm = Zsummary.pres) |>
    drop_na(zsumm) |>
    separate(pair, c("ref_net", "test_net"), sep = "_") |>
    mutate(ref_net = str_remove(ref_net, "ref\\."),
	   test_net = str_remove(test_net, "inColumnsAlsoPresentIn\\.")) |>
    filter(! module %in% c("grey", "gold")) |>
    mutate_at(vars(ref_net, test_net), 
	      ~recode(., "CD40L" = "CD40c", "TLR9" = "TLR9c", "TLR7" = "TLR7c",
		      "BCR" = "BCRc", "BCR-TLR7" = "BCR/TLR7c", "DN2" = "DN2c")) |>
    mutate_at(vars(ref_net, test_net),
	      ~factor(., levels = c("CD40c", "TLR9c", "TLR7c", "BCRc", "BCR/TLR7c", "DN2c"))) |>
    group_by(ref_net, test_net) |>
    arrange(desc(ngenes)) |>
    mutate(module_id = factor(1:n())) |>
    ungroup() |>
    arrange(ref_net, test_net)

z_plot <- 
    ggplot(z_tidy, aes(x = module_id, y = test_net)) +
    geom_tile(aes(fill = zsumm), color = "white") +
    #scale_fill_viridis_c(option = "magma") +
    scale_fill_gradient2(low = "beige",
			mid = "gold",
			high = "tomato3",
			midpoint = 10,
			limits = c(0, 20),
			oob = scales::squish
    ) +
    facet_grid(.~ref_net, scales = "free_x", space = "free_x") +
    theme_minimal() +
    theme(
	  axis.text = element_text(size = 7),
	  axis.title.x = element_text(size = 8),
	  legend.text = element_text(size = 7),
	  legend.title = element_text(size = 8),
	  strip.text = element_text(size = 8),
	  panel.grid.major.x = element_line(color = "gray90", linewidth = .2),
	  panel.grid.major.y = element_blank(),
	  panel.spacing.x = unit(.1, "lines"),
	  legend.box.margin = margin(l = -10),
	  plot.background = element_rect(fill = "white", color = "white")
    ) +
    guides(fill = guide_colorbar(barwidth = .5, barheight = 5)) +
    labs(x = "Module", y = NULL, fill = expression(Z[summ]))

ggsave("./plots/zsummary.png", z_plot, height = 1.75, width = 6.5)
