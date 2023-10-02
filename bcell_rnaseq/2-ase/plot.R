library(tidyverse)
library(qvalue)
library(ggrepel)
library(RColorBrewer)
library(patchwork)
library(ggh4x)

if (!file.exists("plots")) dir.create("plots")

meta <- 
    "./array_spec.tsv" |>
    read_tsv(col_names = c("donor_id", "sample_id", "stim", "mgbid"), 
	     col_types = c(.default = "c")) |>
    unite("id", c("sample_id", "stim"), sep = "_", remove = FALSE) |>
    select(id, donor_id, sample_id, stim)
    
ase_df <- 
    sprintf("./results/%s.asereadcounter.txt", meta$id) |>
    setNames(meta$id) |>
    map_df(read_tsv, .id = "id") |>
    left_join(meta, by = "id") |>
    select(donor_id, sample_id, stim, everything()) |>
    select(-id)

sample_order <- ase_df |> 
    distinct(sample_id, stim) |> 
    count(sample_id, sort = T) |>
    pull(sample_id)

ref_ratios <- ase_df |>
    select(sample_id, stim, var_id = variantID, refCount, totalCount, otherBases, rawDepth) |>
    mutate(sample_id = factor(sample_id, levels = sample_order),
	   stim = recode(stim, "unstday0" = "Day 0"),
	   stim = factor(stim, levels = c("Day 0", "BCR", "TLR7", "DN2")),
	   total = totalCount + otherBases,
	   ref_r = refCount / total) |>
    select(sample_id, stim, var_id, ref_r)

ref_r_plot <- 
    ggplot(ref_ratios, aes(ref_r)) +
    geom_histogram(aes(fill = stim)) +
    scale_x_continuous(breaks = c(0, .5, 1), labels = c("0", "0.5", "1")) +
    scale_fill_manual(values = c("Day 0" = "grey", "BCR" = "cornflowerblue",
				  "TLR7" = "forestgreen", "DN2" = "tomato3")) +
    facet_grid2(sample_id~stim, scales = "free_y", independent = "y") +
    theme_minimal() +
    theme(text = element_text(size = 10),
	  panel.grid = element_blank(),
	  strip.text.y = element_text(angle = 0),
	  strip.background = element_rect(fill = "white", color = "white"),
	  legend.position = "none",
	  axis.text.y = element_blank(),
	  plot.background = element_rect(fill = "white", color = "white")) +
    labs(x = "Reference allele ratio", y = NULL, fill = "Stim:")

ggsave("./plots/ref_r.png", ref_r_plot, height = 7, width = 4)


