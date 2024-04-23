library(tidyverse)
library(aod)
library(lme4)
library(furrr)

if (!file.exists("results_glm")) dir.create("results_glm")

# Functions
run_model <- function(count_data, random = TRUE) {
   
    # format data
    counts_df <- 
	count_data |>
	mutate(ref_allele = map2(ref_count, alt_count, ~c(rep(1, .x), rep(0, .y)))) |>
	select(variant_id, stim, ref_allele) |>
	unnest(cols = ref_allele) |>
	mutate(ref_allele = as.factor(ref_allele),
	       time_num = as.integer(stim) - 1L)

    out_df <- count_data |>
	filter(stim != "Day 0") |>
	select(donor_id, replic, stim, variant_id)

    if (random == FALSE) {
	# Model with no random effects

	## Intercept-only model
	lm0 <- glm(ref_allele ~ 1, data = counts_df, family = binomial)

	## Include time as fixed
	lm1 <- glm(ref_allele ~ 1 + time_num, data = counts_df, family = binomial)

	anova_res <- anova(lm0, lm1, test = "Chisq")
	out_df$p <- anova_res$P[2]
	out_df$beta_time <- coef(lm1)[["time_num"]]

    } else {

	# Model with random effects

	## Null model with intercept and sample as random effect
	lm0 <- glmer(ref_allele ~ 1 + (1 | stim), data = counts_df, family = binomial)

	## Alternative model
	lm1 <- glmer(ref_allele ~ 1 + (1 | stim) + time_num, data = counts_df, family = binomial)

	anova_res <- anova(lm1, lm0)
	out_df$p <- anova_res$P[2]
    }

    list("dat" = out_df, "anova" = anova_res)
}

quietly_run_model_random <- quietly(function(x) run_model(x, random = TRUE))
quietly_run_model_norandom <- quietly(function(x) run_model(x, random = FALSE))

# Data
# Pre-filtered data
ase_data <- 
    read_tsv("./ase_data.tsv") |>
    separate(sample_id, c("donor_id", "replic"), sep = "_") |>
    mutate(replic = LETTERS[as.integer(replic)],
	   replic = factor(replic)) |>
    select(donor_id, replic, stim, variant_id, ref_count, alt_count)

# Separate data per stim
ase_bcr <- ase_data |>
    filter(stim %in% c("Day 0", "BCR")) |>
    mutate(stim = fct_inorder(stim)) |>
    group_by(donor_id, replic, variant_id) |>
    filter(n_distinct(stim) == 2) |>
    ungroup() |>
    select(donor_id, replic, stim, variant_id, ref_count, alt_count)

ase_tlr <- ase_data |>
    filter(stim %in% c("Day 0", "TLR7")) |>
    mutate(stim = fct_inorder(stim)) |>
    group_by(donor_id, replic, variant_id) |>
    filter(n_distinct(stim) == 2) |>
    ungroup() |>
    select(donor_id, replic, stim, variant_id, ref_count, alt_count)

ase_dn2 <- ase_data |>
    filter(stim %in% c("Day 0", "DN2")) |>
    mutate(stim = fct_inorder(stim)) |>
    group_by(donor_id, replic, variant_id) |>
    filter(n_distinct(stim) == 2) |>
    ungroup() |>
    select(donor_id, replic, stim, variant_id, ref_count, alt_count)


## Run model with random effects
#plan(sequential)
#plan(multisession, workers = availableCores())
#
#bcr_res <- 
#    ase_bcr |>
#    group_split(donor_id, replic, variant_id) |>
#    future_map(quietly_run_model_random)
#
#write_rds(bcr_res, "./results_glm/bcr.rds")
#
#tlr_res <- 
#    ase_tlr |>
#    group_split(donor_id, replic, variant_id) |>
#    future_map(quietly_run_model_random)
#
#write_rds(tlr_res, "./results_glm/tlr.rds")
#
#dn2_res <- 
#    ase_dn2 |>
#    group_split(donor_id, replic, variant_id) |>
#    future_map(quietly_run_model_random)
#
#write_rds(dn2_res, "./results_glm/dn2.rds")



# Run model without random effects
plan(sequential)
plan(multisession, workers = availableCores())

bcr_res_nr <- 
    ase_bcr |>
    group_split(donor_id, replic, variant_id) |>
    future_map(quietly_run_model_norandom)

write_rds(bcr_res_nr, "./results_glm/bcr_nr.rds")

tlr_res_nr <- 
    ase_tlr |>
    group_split(donor_id, replic, variant_id) |>
    future_map(quietly_run_model_norandom)

write_rds(tlr_res_nr, "./results_glm/tlr_nr.rds")

dn2_res_nr <- 
    ase_dn2 |>
    group_split(donor_id, replic, variant_id) |>
    future_map(quietly_run_model_norandom)

write_rds(dn2_res_nr, "./results_glm/dn2_nr.rds")

