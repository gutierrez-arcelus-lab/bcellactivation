library(tidyverse)
library(Seurat)
library(glue)

# WGCNA modules
stims <- c("TLR7", "BCR", "DN2")

module_sizes <-
    glue("../bcell_lowinput/wgcna/data/{stims}_modules.tsv") |>
    setNames(stims) |>
    map_dfr(~read_tsv(.) |>
	    count(module) |>
	    arrange(desc(n)) |>
	    filter(module != "grey") |>
	    rowid_to_column("ix") |>
	    select(module, ix) |>
	    mutate(ix = factor(ix)),
	    .id = "stim")

module_df <- 
    glue("../bcell_lowinput/wgcna/data/{stims}_kme.tsv") |>
    setNames(stims) |>
    map_dfr(~read_tsv(.) |>
	    mutate(gene_id = str_remove(gene_id, "\\.\\d+$")) |>
	    pivot_longer(-gene_id, names_to = "module", values_to = "kme") |>
	    filter(kme >= .9, module != "grey") |>
	    group_by(gene_id) |>
	    slice_max(kme) |>
	    ungroup(),
	    .id = "stim")

module_signatures <- 
    module_df |>
    unite("mod", c(stim, module), sep = "_") |>
    select(-kme) |>
    {function(x) split(x, x$mod)}() |>
    map(~pull(., gene_id))

# scRNA data
bcells <- read_rds("../citeseq/data/seuratv4_qced.rds")

# Add module score
bcells_sig <- 
    AddModuleScore(object = bcells, 
		   features = module_signatures,
		   name = "WGCNA_score")
		   
rename_vec <- 
    setNames(glue("WGCNA_score{seq_along(module_signatures)}"), glue("mod_{names(module_signatures)}"))

heat_data <- 
    bcells_sig@meta.data |> 
    as_tibble(rownames = "barcode") |> 
    rename(any_of(rename_vec)) |>
    select(barcode, seurat_clusters, starts_with("mod")) |>
    pivot_longer(starts_with("mod"), names_to = "module", values_to = "score") |>
    group_by(seurat_clusters, module) |>
    summarise(mean_score = mean(score)) |>
    ungroup() |>
    group_by(module) |>
    mutate(z_score = as.numeric(scale(mean_score))) |> 
    ungroup() |>
    mutate(module = str_remove(module, "mod_")) |>
    separate(module, c("stim", "module"), sep = "_") |>
    left_join(module_sizes, join_by(stim, module)) |>
    mutate(stim = recode(stim, "TLR7" = "TLR7c", "BCR" = "BCRc", "DN2" = "DN2c"),
	   stim = factor(stim, levels = c("TLR7c", "BCRc", "DN2c")))

heat <- 
    ggplot(heat_data, aes(x = ix, y = seurat_clusters)) +
    geom_tile(aes(fill = z_score), color = "white") +
    scale_fill_gradient2(low = "Light Sky Blue", mid = "lightyellow", high = "Dark Red",
			 midpoint = 0) +
    facet_grid(.~stim, space = "free_x", scales = "free_x") +
    theme_minimal() +
    theme(
	  axis.text = element_text(size = 7),
	  axis.title = element_text(size = 8),
	  strip.text = element_text(size = 8),
	  legend.text = element_text(size = 7),
	  legend.title = element_text(size = 8),
	  panel.spacing.x = unit(.5, "lines"),
	  legend.box.margin = margin(l = -10),
	  plot.background = element_rect(fill = "white", color = "white")
    ) +
    guides(fill = guide_colorbar(barwidth = .5)) +
    labs(x = "Module", y = "Cluster", fill = "Score")

ggsave("modules_in_scRNAseq.png", heat, height = 2, width = 6.5)
