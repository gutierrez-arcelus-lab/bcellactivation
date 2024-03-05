library(tidyverse)
library(aod)
library(lme4)
library(furrr)

# Functions
run_model <- function(count_data) {
   
    # format data
    counts_df <- 
	count_data |>
	mutate(ref_allele = map2(ref_count, alt_count, ~c(rep(1, .x), rep(0, .y)))) |>
	select(var_id, stim, ref_allele) |>
	unnest(cols = ref_allele) |>
	mutate(ref_allele = as.factor(ref_allele),
	       time_num = as.integer(stim) - 1L)

    # check if intercept is different from zero
    #lm0 <- glm(ref_allele ~ 1, data = counts_df, family = binomial)

    # Null model with intercept and sample as random effect
    lm1 <- glmer(ref_allele ~ 1 + (1 | stim), data = counts_df, family = binomial)

    # alternative model providing stim term
    lm2 <- glmer(ref_allele ~ 1 + (1 | stim) + time_num, data = counts_df, family = binomial)

    list("variant" = unique(counts_df$var_id), "anova" = anova(lm2, lm1))
}

quietly_run_model <- quietly(run_model)

# Data
# Pre-filtered data
ase_data <- 
    read_tsv("./ase_data_20reads.tsv") |>
    separate(sample_id, c("donor_id", "replic"), sep = "\\.") |>
    mutate(replic = LETTERS[as.integer(replic)],
	   replic = factor(replic)) |>
    select(donor_id, replic, stim, var_id, ref_count, alt_count)

# Separate data per stim
ase_bcr <- ase_data |>
    filter(stim %in% c("Day 0", "BCR")) |>
    mutate(stim = fct_inorder(stim)) |>
    group_by(donor_id, replic, var_id) |>
    filter(n_distinct(stim) == 2) |>
    ungroup() |>
    select(donor_id, replic, stim, var_id, ref_count, alt_count)

ase_tlr <- ase_data |>
    filter(stim %in% c("Day 0", "TLR7")) |>
    mutate(stim = fct_inorder(stim)) |>
    group_by(donor_id, replic, var_id) |>
    filter(n_distinct(stim) == 2) |>
    ungroup() |>
    select(donor_id, replic, stim, var_id, ref_count, alt_count)

ase_dn2 <- ase_data |>
    filter(stim %in% c("Day 0", "DN2")) |>
    mutate(stim = fct_inorder(stim)) |>
    group_by(donor_id, replic, var_id) |>
    filter(n_distinct(stim) == 2) |>
    ungroup() |>
    select(donor_id, replic, stim, var_id, ref_count, alt_count)

plan(multisession, workers = availableCores())

bcr_res <- 
    ase_bcr |>
    group_split(donor_id, replic, var_id) |>
    future_map(quietly_run_model)

write_rds(bcr_res, "./results_glm/bcr.rds")

tlr_res <- 
    ase_tlr |>
    group_split(donor_id, replic, var_id) |>
    future_map(quietly_run_model)

write_rds(tlr_res, "./results_glm/tlr.rds")

dn2_res <- 
    ase_dn2 |>
    group_split(donor_id, replic, var_id) |>
    future_map(quietly_run_model)

write_rds(dn2_res, "./results_glm/dn2.rds")



#var_ids <- bcr_res |>
#    map("result") |>
#    map_chr("variant")
#
#pvals <- bcr_res |>
#    map("result") |>
#    map("anova") |> 
#    map_dbl(function(x) x$P[2])
#
#tibble(var_id = var_ids, p = pvals) |>
#    arrange(p)
#

