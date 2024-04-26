library(tidyverse)
library(aod)
library(lme4)
library(furrr)

if (!file.exists("results_glm")) dir.create("results_glm")

# Functions
run_model <- function(count_data, random = FALSE) {
   
    # format data
    counts_df <- 
	count_data |>
	mutate(ref_allele = map2(ref_count, alt_count, ~c(rep(1, .x), rep(0, .y)))) |>
	select(variant_id, stim, ref_allele) |>
	unnest(cols = ref_allele) |>
	mutate(ref_allele = as.factor(ref_allele),
	       time_num = as.integer(stim) - 1L)

    out_df <- 
	count_data |>
	group_by(donor_id, replic, variant_id) |>
	summarise(stim = paste(stim, collapse = "-")) |>
	ungroup() |>
	select(donor_id, replic, stim, variant_id)

    
    if (isTRUE(random)) {

	# Model with random effects

	## Null model with intercept and sample as random effect
	lm0 <- glmer(ref_allele ~ 1 + (1 | stim), data = counts_df, family = binomial)

	## Alternative model
	lm1 <- glmer(ref_allele ~ 1 + (1 | stim) + time_num, data = counts_df, family = binomial)

	anova_res <- anova(lm1, lm0)
	out_df$p <- anova_res$P[2]
    }

    if (isFALSE(random)) {
	# Model with no random effects

	## Intercept-only model
	lm0 <- glm(ref_allele ~ 1, data = counts_df, family = binomial)

	## Include time as fixed
	lm1 <- glm(ref_allele ~ 1 + time_num, data = counts_df, family = binomial)

	anova_res <- anova(lm0, lm1, test = "Chisq")
	out_df$p <- anova_res$P[2]
	out_df$beta_time <- coef(lm1)[["time_num"]]

    }

    list("dat" = out_df, "anova" = anova_res)
}

quietly_run_model <- quietly(function(x) run_model(x, random = FALSE))

filter_data <- function(stim_a, stim_b) {

    ase_data |>
	filter(stim %in% c(stim_a, stim_b)) |>
	mutate(stim = factor(stim, levels = c(stim_a, stim_b))) |>
	group_by(donor_id, replic, variant_id) |>
	filter(n_distinct(stim) == 2) |>
	ungroup() |>
	select(donor_id, replic, stim, variant_id, ref_count, alt_count)
}

# Data
# Pre-filtered data
stim_levels <- c("Day 0", "BCR", "TLR7", "DN2") 

ase_data <- 
    read_tsv("./ase_data.tsv") |>
    separate(sample_id, c("donor_id", "replic"), sep = "_") |>
    mutate(replic = LETTERS[as.integer(replic)],
	   replic = factor(replic)) |>
    select(donor_id, replic, stim, variant_id, ref_count, alt_count)

		       
# Run model without random effects for all stim pairs
stim_pairs <- 
    ase_data |>
    distinct(stim_1 = stim) |>
    expand_grid(stim_2 = stim_1) |>
    mutate_at(vars(stim_1, stim_2), ~factor(., levels = stim_levels)) |>
    filter(as.numeric(stim_1) < as.numeric(stim_2)) |>
    mutate_at(vars(stim_1, stim_2), as.character)

plan(sequential)
plan(multisession, workers = availableCores())

res_df <- 
    stim_pairs |>
    filter(stim_1 == "Day 0") |>
    mutate(data = map2(stim_1, stim_2, 
		       ~filter_data(.x, .y) |>
		       group_split(donor_id, replic, variant_id) |>
		       future_map(quietly_run_model)))

res_df$data |>
    flatten() |>
    write_rds("./results_glm/glm_res.rds")

res_df$data |>
    flatten() |>
    map("result") |>
    map_df("dat") |>
    write_tsv("./results_glm/glm_res_df.tsv")



#bcr_res_nr <- 
#    ase_bcr |>
#    group_split(donor_id, replic, variant_id) |>
#    future_map(quietly_run_model_norandom)
#
#write_rds(bcr_res_nr, "./results_glm/bcr_nr.rds")
#
#tlr_res_nr <- 
#    ase_tlr |>
#    group_split(donor_id, replic, variant_id) |>
#    future_map(quietly_run_model_norandom)
#
#write_rds(tlr_res_nr, "./results_glm/tlr_nr.rds")
#
#dn2_res_nr <- 
#    ase_dn2 |>
#    group_split(donor_id, replic, variant_id) |>
#    future_map(quietly_run_model_norandom)
#
#write_rds(dn2_res_nr, "./results_glm/dn2_nr.rds")
#
