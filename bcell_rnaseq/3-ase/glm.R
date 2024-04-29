library(tidyverse)
library(aod)
library(lme4)
library(furrr)

if (!file.exists("results_glm")) dir.create("results_glm")

# Functions
filter_data <- function(stim_a, stim_b) {

    ase_data |>
	filter(stim %in% c(stim_a, stim_b)) |>
	mutate(stim = factor(stim, levels = c(stim_a, stim_b))) |>
	group_by(donor_id, replic, variant_id) |>
	filter(n_distinct(stim) == 2) |>
	ungroup() |>
	select(donor_id, replic, stim, variant_id, ref_count, alt_count)
}

run_model <- function(count_data) {

    counts_df <-
	count_data |>
	mutate(total = ref_count + alt_count,
	       ref_r = ref_count/total) |>
	select(stim, ref_r, total)
    
    out_df <- 
	count_data |>
	group_by(donor_id, replic, variant_id) |>
	summarise(stim = paste(sort(stim), collapse = "-")) |>
	ungroup() |>
	select(donor_id, replic, stim, variant_id)
    
    lm1 <- glm(ref_r ~ stim, data = counts_df, family = binomial, weights = total) 
    anova_res <- anova(lm1, test = "Chisq")
    out_df$p <- anova_res$P[2]
    out_df$beta_stim <- coef(lm1)[[2]]

    list("dat" = out_df, "anova" = anova_res)
}

quietly_run_model <- quietly(function(x) run_model(x))

# Data
# Pre-filtered data
stim_levels <- c("Day 0", "BCR", "TLR7", "DN2") 

ase_data <- 
    read_tsv("./ase_data.tsv") |>
    separate(sample_id, c("donor_id", "replic"), sep = "_") |>
    mutate(replic = LETTERS[as.integer(replic)],
	   replic = factor(replic),
	   stim = factor(stim, levels = stim_levels)) |>
    select(donor_id, replic, stim, variant_id, ref_count, alt_count) |>
    arrange(donor_id, replic, stim, variant_id)
		       
# Run model without random effects for all stim pairs
stim_pairs <- 
    ase_data |>
    distinct(stim_1 = stim) |>
    expand_grid(stim_2 = stim_1) |>
    filter(as.numeric(stim_1) < as.numeric(stim_2)) |>
    mutate_at(vars(stim_1, stim_2), as.character)

plan(multisession, workers = availableCores())

res_df <- 
    stim_pairs |>
    mutate(data = map2(stim_1, stim_2, 
		       ~filter_data(.x, .y) |>
		       group_split(donor_id, replic, variant_id) |>
		       future_map(quietly_run_model)))

# Save all output
res_df$data |>
    flatten() |>
    write_rds("./results_glm/glm_res.rds")

# Save coefficients
res_df$data |>
    flatten() |>
    map("result") |>
    map_df("dat") |>
    write_tsv("./results_glm/glm_res_df.tsv")

