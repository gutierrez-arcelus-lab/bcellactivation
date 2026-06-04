# ==============================================================================
# Description:  Evaluates the robustness and reproducibility of gene co-expression 
#               modules across different stimulation conditions. Uses permutation 
#               testing to calculate the Z-summary statistic for all-vs-all 
#               network comparisons.
# Input:        1. *_counts.rds (Expression matrices for all 6 conditions)
#               2. *_modules.tsv (Module assignments for all 6 conditions)
# Output:       1. module_preservation.rds (Raw WGCNA permutation results)
#               2. fig_s3.png (Heatmap of Z-summary statistics)
# ==============================================================================

library(tidyverse)
library(WGCNA)

# ------------------------------------------------------------------------------
# 1. Data Preparation
# ------------------------------------------------------------------------------
# WGCNA requires a specific nested list structure for multi-network analysis.
# The top level is the condition name, containing a list with the 'data' matrix.
multi_expr <-
    list(
	 "CD40L"    = list(data = read_rds("../01_rnaseq_lowinput/3_wgcna/data/CD40L_counts.rds")),
	 "TLR9"     = list(data = read_rds("../01_rnaseq_lowinput/3_wgcna/data/TLR9_counts.rds")),
	 "TLR7"     = list(data = read_rds("../01_rnaseq_lowinput/3_wgcna/data/TLR7_counts.rds")),
	 "BCR"      = list(data = read_rds("../01_rnaseq_lowinput/3_wgcna/data/BCR_counts.rds")),
	 "BCR-TLR7" = list(data = read_rds("../01_rnaseq_lowinput/3_wgcna/data/BCR-TLR7_counts.rds")),
	 "DN2"      = list(data = read_rds("../01_rnaseq_lowinput/3_wgcna/data/DN2_counts.rds"))
    )

# Create a matching list of named vectors containing the module color assignments 
# for each gene in each condition.
multi_color <- 
    list(
	 "CD40L"    = read_tsv("../01_rnaseq_lowinput/3_wgcna/data/CD40L_modules.tsv") |> deframe(),
	 "TLR9"     = read_tsv("../01_rnaseq_lowinput/3_wgcna/data/TLR9_modules.tsv") |> deframe(),
	 "TLR7"     = read_tsv("../01_rnaseq_lowinput/3_wgcna/data/TLR7_modules.tsv") |> deframe(),
	 "BCR"      = read_tsv("../01_rnaseq_lowinput/3_wgcna/data/BCR_modules.tsv") |> deframe(),
	 "BCR-TLR7" = read_tsv("../01_rnaseq_lowinput/3_wgcna/data/BCR-TLR7_modules.tsv") |> deframe(),
	 "DN2"      = read_tsv("../01_rnaseq_lowinput/3_wgcna/data/DN2_modules.tsv") |> deframe()
    )

# ------------------------------------------------------------------------------
# 2. Module Preservation Permutation Testing
# ------------------------------------------------------------------------------
# Calculate preservation statistics across all networks.
# - referenceNetworks = 1:6 instructs WGCNA to do an all-vs-all comparison 
#   (every condition acts as a reference against all other test conditions).
# - nPermutations = 100 randomly shuffles the data 100 times to build a null 
#   distribution for calculating the Z-summary significance.
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

# ------------------------------------------------------------------------------
# 3. Parse and Reshape Results
# ------------------------------------------------------------------------------
mp <- read_rds("data/module_preservation.rds")


# Extract the Z-statistics from the heavily nested WGCNA output object.
# General interpretation of Zsummary.pres:
# > 10 : Highly preserved (the module exists in both conditions)
# 2 to 10 : Weakly to moderately preserved
# < 2 : Not preserved (the module is unique to the reference condition)
z_tidy <- 
    mp$preservation$Z |>
    list_flatten() |>
    map_dfr(~as_tibble(., rownames = "module"), .id = "pair") |>
    # Select only the composite Z-summary metric and the module size
    select(pair, module, ngenes = moduleSize, zsumm = Zsummary.pres) |>
    drop_na(zsumm) |>
    # Clean up WGCNA's automated 'pair' naming convention 
    # (e.g., "ref.CD40L_inColumnsAlsoPresentIn.TLR9")
    separate(pair, c("ref_net", "test_net"), sep = "_") |>
    mutate(ref_net = str_remove(ref_net, "ref\\."),
	   test_net = str_remove(test_net, "inColumnsAlsoPresentIn\\.")) |>
    # Remove unassigned genes (grey) and the random sampling module (gold)
    filter(! module %in% c("grey", "gold")) |>
    # Recode stimulation names to match the final manuscript figure labels
    mutate_at(vars(ref_net, test_net), 
	      ~recode(., 
		      "CD40L" = "CD40c", 
		      "TLR9" = "TLR9c", 
		      "TLR7" = "TLR7c",
		      "BCR" = "BCRc", 
		      "BCR-TLR7" = "BCR/TLR7c", 
		      "DN2" = "DN2c")) |>
    mutate_at(vars(ref_net, test_net),
	      ~factor(., levels = c("CD40c", "TLR9c", "TLR7c", "BCRc", "BCR/TLR7c", "DN2c"))) |>
    # Assign an arbitrary numeric ID to modules based on size for plotting
    group_by(ref_net, test_net) |>
    arrange(desc(ngenes)) |>
    mutate(module_id = factor(1:n())) |>
    ungroup() |>
    arrange(ref_net, test_net)

# ------------------------------------------------------------------------------
# 4. Plot Z-Summary Heatmap (Figure S3)
# ------------------------------------------------------------------------------
z_plot <- 
    ggplot(z_tidy, aes(x = module_id, y = test_net)) +
    geom_tile(aes(fill = zsumm), color = "white") +
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

ggsave("./sfigs/sfig3_zsummary.png", z_plot, height = 1.75, width = 6.5)
