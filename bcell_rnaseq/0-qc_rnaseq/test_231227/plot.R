library(tidyverse)
library(glue)
library(ggrepel)

array_spec <- 
    read_tsv("./array_spec_mbv.tsv", col_names = c("prefix", "batch"))

mbv_df <- 
    array_spec |>
    unite("prefix", c("prefix", "batch"), sep = "_") |>
    mutate(f = glue("results/{prefix}.txt")) |>
    deframe() |>
    map_dfr(read_table, .id = "prefix") |> 
    separate(prefix, c("rnaseq_donor_id", "rnaseq_rep_id", "rnaseq_stim", "vcf_batch"), sep = "_") |>
    unite("rnaseq_sample_id", c(rnaseq_donor_id, rnaseq_rep_id), sep = "_", remove = FALSE) |>
    mutate(vcf_donor_id = str_extract(SampleID, ".+-(\\d+)$", group = 1)) |>
    select(starts_with("rnaseq"), 
	   vcf_batch, 
	   vcf_donor_id, 
	   vcf_sample_id = SampleID, 
	   n_het_covered, 
	   n_hom_covered, 
	   perc_het_consistent,
	   perc_hom_consistent)

label_df <- 
    mbv_df |>
    group_by(rnaseq_sample_id, rnaseq_stim) |>
    slice_max(perc_het_consistent + perc_hom_consistent) |>
    ungroup() |>
    filter(rnaseq_donor_id != vcf_donor_id)
    

test_p <- 
    ggplot(mbv_df, 
	   aes(x = perc_het_consistent, y = perc_hom_consistent)) +
    geom_point(aes(color = rnaseq_donor_id == vcf_donor_id)) +
    scale_x_continuous(limits = c(0, 1),
		       breaks = c(0, .5, 1),
		       labels = c("0", "0.5", "1")) +
    scale_y_continuous(limits = c(0, 1),
		       breaks = c(0, .5, 1),
		       labels = c("0", "0.5", "1")) +
    scale_color_manual(values = c("TRUE" = "Lime Green", "FALSE" = "tomato3")) +
    geom_text_repel(data = label_df, 
		    aes(x = perc_het_consistent, 
			y = perc_hom_consistent, 
			label = vcf_donor_id),
		    direction = "y",
		    hjust = "inward",
		    nudge_y = -.5,
		    segment.size = .2,
		    size = 3) +
    facet_wrap(rnaseq_sample_id ~ rnaseq_stim) +
    theme_bw() +
    theme(panel.grid = element_blank(),
	  legend.margin = margin(0, 0, 0, 0),
	  legend.position = "top") +
    labs(x = "Proportion consistent Hets", 
	 y = "Proportion consistent Homs",
	 color = "Match:") +
    guides(color = guide_legend(override.aes = list(size = 2)))

ggsave("./mbv.png", test_p, width = 6, height = 7)
