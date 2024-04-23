library(tidyverse)
library(glue)
library(ggrepel)

array_spec <- 
    read_tsv("./array_spec_mbv.txt", col_names = c("prefix", "batch"))

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
	   perc_hom_consistent) |>
    mutate(rnaseq_stim = recode(rnaseq_stim, "unstday0" = "Day 0"),
	   rnaseq_stim = factor(rnaseq_stim, levels = c("Day 0", "BCR", "TLR7", "DN2")))

label_df <- 
    mbv_df |>
    group_by(rnaseq_sample_id, rnaseq_stim) |>
    slice_max(perc_het_consistent + perc_hom_consistent) |>
    ungroup() |>
    filter(rnaseq_donor_id != vcf_donor_id)

mbv_plot <- 
    ggplot(data = mbv_df,
	   aes(x = perc_het_consistent, 
	       y = perc_hom_consistent,
	       color = rnaseq_donor_id == vcf_donor_id)) +
    geom_point(data = filter(mbv_df, rnaseq_donor_id != vcf_donor_id),
	       size = 2) +
    geom_point(data = filter(mbv_df, rnaseq_donor_id == vcf_donor_id),
	       size = 2) +
    scale_x_continuous(limits = c(0, 1),
		       breaks = c(0, 1),
		       labels = c("0", "1")) +
    scale_y_continuous(limits = c(0, 1),
		       breaks = c(0, 1),
		       labels = c("0", "1")) +
    scale_color_manual(values = c("TRUE" = "Yellow Green", "FALSE" = "Dark Red")) +
    geom_text_repel(data = label_df, 
		    aes(x = perc_het_consistent, 
			y = perc_hom_consistent, 
			label = vcf_donor_id),
		    direction = "y",
		    hjust = "inward",
		    nudge_y = -.5,
		    segment.size = .2,
		    size = 3) +
    facet_grid(rnaseq_stim ~ sub("_", "\n", rnaseq_sample_id)) +
    theme_bw() +
    theme(axis.text = element_text(size = 10),
	  axis.title = element_text(size = 12),
	  panel.grid = element_blank(),
	  panel.spacing = unit(0.2, "lines"),
	  strip.text.x = element_text(size = 10),
	  strip.text.y = element_text(angle = 0, size = 11),
	  strip.background = element_rect(fill = "white", color = "white")) +
    labs(x = "Proportion consistent Heterozygous sites", 
	 y = "Proportion consistent Homs") +
    guides(color = "none")

ggsave("./mbv.png", mbv_plot, width = 15, height = 3.15)


## Save matching results
#mbv_df |>
#    group_by(vcf_donor_id) |>
#    mutate(n_batches = n_distinct(vcf_batch)) |>
#    filter(n_batches == 1 | vcf_batch != "0410") |>
#    group_by(rnaseq_sample_id, rnaseq_stim) |>
#    slice_max(perc_het_consistent + perc_hom_consistent) |>
#    ungroup() |>
#    select(rnaseq_sample_id, rnaseq_stim, vcf_donor_id, vcf_sample_id, vcf_batch) |>
#    write_tsv("matching_results.tsv")
#
