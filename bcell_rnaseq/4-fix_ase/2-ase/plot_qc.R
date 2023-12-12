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
    select(-id) |>
    mutate(stim = recode(stim, "unstday0" = "Day 0"),
	   stim = factor(stim, levels = c("Day 0", "BCR", "TLR7", "DN2")))

sample_order <- ase_df |> 
    distinct(sample_id, stim) |> 
    count(sample_id, sort = TRUE) |>
    pull(sample_id)

ase_clean_df <- ase_df |>
    filter(totalCount >= 10,
	   otherBases/refCount < .1 & otherBases/altCount < .1,
	   (otherBases/(refCount + altCount)) < 0.05) |>
    select(sample_id, stim, var_id = variantID, refCount, altCount) |>
    mutate(sample_id = factor(sample_id, levels = sample_order))

write_tsv(ase_clean_df, "./ase_data_clean.tsv")


ref_ratios <- ase_clean_df |>
    mutate(ref_r = refCount / (refCount + altCount)) |>
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

# Fraction of both alleles seen
both_seen_df <- ase_df |>
    filter(totalCount >= 10) |>
    mutate(both_seen = refCount >= 1 & altCount >= 1) |>
    mutate(sample_id = factor(sample_id, levels = sample_order),
	   stim = recode(stim, "unstday0" = "Day 0"),
	   stim = factor(stim, levels = c("Day 0", "BCR", "TLR7", "DN2"))) |>
    group_by(donor_id, sample_id, stim) |>
    summarise(p_both_seen = mean(both_seen)) |>
    ungroup() 

both_seen_plot <- 
    ggplot(both_seen_df |> mutate(sample_id = fct_rev(sample_id)), 
	   aes(x = p_both_seen, y = sample_id)) +
	geom_col(aes(fill = stim)) +
	scale_x_continuous(limits = c(0, 1),
			   breaks = c(0, .5, 1),
			   labels = c("0", "0.5", "1")) +
	scale_fill_manual(values = c("Day 0" = "grey", "BCR" = "cornflowerblue",
				     "TLR7" = "forestgreen", "DN2" = "tomato3")) +
	facet_wrap(~stim, nrow = 1) +
	theme_minimal() +
	theme(panel.grid.minor.x = element_blank(),
	      legend.position = "none",
	      plot.background = element_rect(fill = "white", color = "white")) +
	labs(x = "Proportion of both alleles seen", y = "Sample ID")

ggsave("./plots/both_seen.png", both_seen_plot)

## Annotate genes
#annotations <- 
#    file.path("/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens",
#	      "gencode.v39.primary_assembly.annotation.gtf") |>
#    read_tsv(comment = "#", col_names = FALSE, col_types = "ccciicccc")
#
#gene_df <- annotations |>
#    filter(X3 == "gene", X1 %in% paste0("chr", c(1:22, "X"))) |>
#    mutate(gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
#	   gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+")) |>
#    select(chr = X1, start = X4, end = X5, gene_id, gene_name)
#
#ase_res_annot <-
#    left_join(ase_df, gene_df, join_by(contig == chr, between(position, start, end))) |>
#    select(-start, -end)
#
#
#

# compare technical reps
ase_clean_reps <- 
    ase_clean_df |>
    separate(sample_id, c("donor_id", "replic"), sep = "\\.") |>
    group_by(donor_id) |>
    filter(any(replic %in% c(2,3))) |>
    ungroup() |>
    mutate(ref_r = refCount/totalCount,
	   replic = paste0("rep", replic)) |>
    select(donor_id, replic, stim, var_id, ref_r)

total_counts <- ase_df |>
    filter(donor_id == unique(ase_clean_reps$donor_id)[2]) |>
    select(sample_id, stim, var_id = variantID, totalCount) |>
    separate(sample_id, c("donor_id", "replic"), sep = "\\.") |>
    mutate(replic = paste0("total_", replic)) |>
    pivot_wider(names_from = replic, values_from = totalCount) |>
    mutate_at(vars(starts_with("total_")), ~replace_na(., 0)) |>
    select(donor_id, stim, var_id, starts_with("total_"))

temp_df <- ase_clean_reps |>
    filter(donor_id == unique(donor_id)[2]) |>
    pivot_wider(names_from = replic, values_from = ref_r) |>
    select(donor_id, stim, var_id, rep1, rep2) |>
    filter(!(is.na(rep1) & is.na(rep2))) |>
    mutate_at(vars(starts_with("rep")), ~replace_na(., -0.1)) |>
    left_join(total_counts, join_by(donor_id, stim, var_id))

test_p <- 
    ggplot(temp_df, aes(rep1, rep2)) +
    geom_point(aes(color = log10(total_1)), size = .25) +
    scale_x_continuous(breaks = c(-.1, 0, .5, 1),
		       labels = c("NA", "0", "0.5", "1")) +
    scale_y_continuous(breaks = c(-.1, 0, .5, 1),
		       labels = c("NA", "0", "0.5", "1")) +
    scale_color_continuous(low = "grey96", high = "black",
			   guide = guide_colorbar(barheight = .5)) +
    facet_grid(donor_id~stim) +
    theme_bw() +
    theme(panel.grid.minor.x = element_blank(),
	  panel.grid.minor.y = element_blank(),
	  strip.text.y = element_text(angle = 0),
	  legend.position = "top") +
    labs(x = "REF allele ratio in rep #1",
	 y = "REF allele ratio in rep #2",
	 color = "Total counts in rep #1 (log10):")

library(ggpointdensity)

test_p2 <- 
    ggplot(temp_df, aes(rep1, rep2)) +
    geom_pointdensity(size = .25) +
    scale_x_continuous(breaks = c(-.1, 0, .5, 1),
		       labels = c("NA", "0", "0.5", "1")) +
    scale_y_continuous(breaks = c(-.1, 0, .5, 1),
		       labels = c("NA", "0", "0.5", "1")) +
    scale_color_continuous(low = "beige", high = "red",
			   guide = guide_colorbar(barheight = .5, barwidth = 10)) +
    facet_grid(donor_id~stim) +
    theme_bw() +
    theme(panel.grid.minor.x = element_blank(),
	  panel.grid.minor.y = element_blank(),
	  strip.text.y = element_text(angle = 0),
	  legend.position = "top") +
    labs(x = "REF allele ratio in rep #1",
	 y = "REF allele ratio in rep #2",
	 color = "Density:")


## color by significance at binomial test
library(qvalue)

ase_res <- 
    ase_clean_df |>
    mutate(p_value = map2_dbl(refCount, totalCount, 			      
			      ~binom.test(.x, .y, p = .5, alternative = "two.sided")$p.value),
	   q_value = qvalue(p_value)$qvalues) 

temp_pvals <- 
    ase_res |>
    separate(sample_id, c("donor_id", "replic"), sep = "\\.") |>
    filter(donor_id == unique(temp_df$donor_id), replic %in% c(1, 2)) |>
    mutate(replic = paste0("q_", replic)) |>
    select(donor_id, replic, stim, var_id, q_value) |>
    pivot_wider(names_from = replic, values_from = q_value) |>
    mutate_at(vars(q_1, q_2), ~replace_na(., 1)) |>
    mutate(sig = case_when(q_1 <= 0.05 & q_2 > 0.05 ~ "rep #1",
			   q_2 <= 0.05 & q_1 > 0.05 ~ "rep #2",
			   q_1 <= 0.05 & q_2 <= 0.05 ~ "both",
			   TRUE ~ "none")) |>
    select(donor_id, stim, var_id, sig)

test_p3 <- temp_df |>
    left_join(temp_pvals, join_by(donor_id, stim, var_id)) |>
    ggplot(aes(rep1, rep2)) +
    geom_point(aes(color = sig), size = .25) +
    scale_x_continuous(breaks = c(-.1, 0, .5, 1),
		       labels = c("NA", "0", "0.5", "1")) +
    scale_y_continuous(breaks = c(-.1, 0, .5, 1),
		       labels = c("NA", "0", "0.5", "1")) +
    scale_color_manual(values = c("none" = "grey90", "rep #1" = "blue", "rep #2" = "red", both = "green3")) +
    facet_grid(donor_id~stim) +
    theme_bw() +
    theme(panel.grid.minor.x = element_blank(),
	  panel.grid.minor.y = element_blank(),
	  strip.text.y = element_text(angle = 0),
	  legend.position = "top") +
    labs(x = "REF allele ratio in rep #1",
	 y = "REF allele ratio in rep #2",
	 color = "Significant in:") +
    guides(color = guide_legend(override.aes = list(size = 2)))

ggsave("./plots/reps.png", 
       test_p + test_p2 + test_p3 + plot_layout(ncol = 1), 
       height = 8, width = 8)

temp_df |>
    left_join(temp_pvals, join_by(donor_id, stim, var_id)) |>
    count(sig, rep1 == -0.1, rep2 == -0.1)


# For seminar
average_imb <- ref_ratios |>
    group_by(sample_id, stim) |>
    summarise(mean_imb = mean(ref_r)) |>
    ungroup()

ref_r_plot_wide <- 
    ggplot(ref_ratios, aes(ref_r)) +
    geom_histogram(aes(fill = stim)) +
    scale_x_continuous(breaks = c(0, .5, 1), labels = c("0", "0.5", "1")) +
    scale_fill_manual(values = c("Day 0" = "grey", "BCR" = "cornflowerblue",
				  "TLR7" = "forestgreen", "DN2" = "tomato3")) +
    geom_vline(data = average_imb, 
	       aes(xintercept = mean_imb), 
	       linewidth = .5) +
    facet_grid2(stim~sample_id, scales = "free_y", independent = "y") +
    theme_minimal() +
    theme(text = element_text(size = 10),
	  panel.grid = element_blank(),
	  strip.text.x = element_text(size = 6),
	  strip.text.y = element_text(angle = 0, size = 12, face = "bold"),
	  strip.background = element_rect(fill = "white", color = "white"),
	  legend.position = "none",
	  axis.text.y = element_blank(),
	  plot.background = element_rect(fill = "white", color = "white")) +
    labs(x = "Reference allele ratio", y = NULL, fill = "Stim:")

ggsave("./plots/ref_r_wide.png", ref_r_plot_wide, height = 3, width = 12)

