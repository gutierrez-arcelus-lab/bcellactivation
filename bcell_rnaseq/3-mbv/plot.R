library(tidyverse)
library(ggrepel)

if (!file.exists("plots")) dir.create("plots")

# list of candidate individuals sent to MGB for recruitment
recruit <- "../../mgb_biobank/candidate_individuals_ids.txt" |>
    read_lines()

meta <- read_tsv("./array_spec_mbv.txt", col_names = c("prefix", "batch"))

prefix <- paste(meta$prefix, meta$batch, sep = "_")

files <- 
    sprintf("results/%s.txt", prefix) |>
    setNames(prefix)

res <- 
    files |> 
    map_dfr(~read_delim(., delim = " ") |> 
	    separate(SampleID, c("vcf_prefix", "vcf_donor_id"), sep = "-") |>
	    filter(vcf_donor_id %in% recruit),
	    .id = "id") |>
    extract(id, c("bam_id", "bam_stim", "vcf_batch"), "(.+)_(.+)_(.+)") |>
    select(bam_id, bam_stim, vcf_donor_id, vcf_prefix, vcf_batch, 
	   n_het = n_het_covered, n_hom = n_hom_covered, 
	   consist_het = perc_het_consistent, consist_hom = perc_hom_consistent, 
	   n_ase = n_het_in_ase) |>
    separate(bam_id, c("bam_donor_id", "bam_rep"), sep = "\\.")

matches_df <-
    res |>
    filter(bam_donor_id == vcf_donor_id) |>
    select(bam_donor_id, bam_rep, bam_stim, vcf_donor_id, consist_het, consist_hom)

recruit_df <- 
    filter(res, vcf_donor_id %in% recruit) |>
    select(bam_donor_id, bam_rep, bam_stim, vcf_donor_id, consist_het, consist_hom)

plot_df <- 
    bind_rows(matches_df, recruit_df) |>
    distinct() |>
    mutate(bam_stim = recode(bam_stim, "unstday0" = "Day 0"),
	   bam_stim = factor(bam_stim, levels = c("Day 0", "BCR", "TLR7", "DN2"))) |>
    unite("sample_id", c(bam_donor_id, bam_rep), sep = ".", remove = FALSE) |>
    group_by(sample_id, bam_stim) |>
    nest() |>
    ungroup() |>
    add_count(sample_id) |>
    arrange(desc(n), sample_id) |>
    select(-n) |>
    mutate(sample_id = fct_inorder(sample_id)) |>
    unnest(cols = data) |>
    mutate(id_match =  bam_donor_id == vcf_donor_id)

mis_df <- plot_df |>
    group_by(sample_id, bam_stim) |>
    slice_max(consist_het + consist_hom) |>
    ungroup() |>
    filter(bam_donor_id != vcf_donor_id)
  
test_p <- 
    ggplot(plot_df, aes(consist_het, consist_hom)) +
    geom_point(data = filter(plot_df, id_match == FALSE),
	       size = 3, shape = 3,
	       aes(color = id_match)) +  
    geom_point(data = filter(plot_df, id_match == TRUE),
	       size = 3, shape = 3,
	       aes(color = id_match)) +  
    geom_text_repel(data = mis_df,
		    aes(consist_het, consist_hom, label = vcf_donor_id),
		    size = 2.5,
		    direction = "y",
		    nudge_y = -.6,
		    segment.size = .25,
		    min.segment.length = 0) +
    scale_x_continuous(limits = c(0, 1.1), breaks = c(0, 1)) +
    scale_y_continuous(limits = c(0, 1.1), breaks = c(0, 1)) +
    scale_color_manual(values = c("TRUE" = "forestgreen", "FALSE" = "tomato3")) +
    facet_grid(sample_id~bam_stim) +
    theme_bw() +
    theme(text = element_text(size = 14),
	  legend.position = "none",
	  panel.grid.minor.x = element_blank(),
	  panel.grid.minor.y = element_blank(),
	  strip.text.y = element_text(angle = 0),
	  strip.background = element_rect(fill = "white", color = "white"),
	  plot.background = element_rect(fill = "white", color = "white")) +
    labs(x = "Fraction concordant heterozygous sites",
	 y = "Fraction concordant homozygous sites")

ggsave("./plots/mbv.png", test_p, width = 5, height = 12)

out <-
    plot_df |>
    group_by(sample_id, bam_stim) |>
    slice_max(consist_het + consist_hom) |>
    ungroup() |>
    mutate(bam_stim = fct_recode(bam_stim, "unstday0" = "Day 0")) |>
    select(bam_donor_id, sample_id, bam_stim, vcf_donor_id) |>
    left_join(distinct(res, vcf_donor_id, vcf_prefix, vcf_batch)) |>
    filter(vcf_batch != "0410") |>
    unite("vcf_id", c(vcf_prefix, vcf_donor_id), sep = "-", remove = FALSE) |>
    select(bam_donor_id, sample_id, bam_stim, vcf_id, vcf_donor_id, vcf_batch)

write_tsv(out, "./match_summary.tsv")

out |> 
    mutate(recruit = vcf_donor_id %in% recruit) |>
    print(n = Inf)

