library(tidyverse)
library(aod)
library(lme4)
library(furrr)

# Functions
run_model <- function(count_data, random = TRUE) {
   
    # format data
    counts_df <- 
	count_data |>
	mutate(ref_allele = map2(ref_count, alt_count, ~c(rep(1, .x), rep(0, .y)))) |>
	select(var_id, stim, ref_allele) |>
	unnest(cols = ref_allele) |>
	mutate(ref_allele = as.factor(ref_allele),
	       time_num = as.integer(stim) - 1L)

    out_df <- count_data |>
	filter(stim != "Day 0") |>
	select(donor_id, replic, stim, var_id)

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


# Run model with random effects
plan(sequential)
plan(multisession, workers = availableCores())

bcr_res <- 
    ase_bcr |>
    group_split(donor_id, replic, var_id) |>
    future_map(quietly_run_model_random)

write_rds(bcr_res, "./results_glm/bcr.rds")

tlr_res <- 
    ase_tlr |>
    group_split(donor_id, replic, var_id) |>
    future_map(quietly_run_model_random)

write_rds(tlr_res, "./results_glm/tlr.rds")

dn2_res <- 
    ase_dn2 |>
    group_split(donor_id, replic, var_id) |>
    future_map(quietly_run_model_random)

write_rds(dn2_res, "./results_glm/dn2.rds")



# Run model without random effects
plan(sequential)
plan(multisession, workers = availableCores())
bcr_res_nr <- 
    ase_bcr |>
    group_split(donor_id, replic, var_id) |>
    future_map(quietly_run_model_norandom)

write_rds(bcr_res_nr, "./results_glm/bcr_nr.rds")

tlr_res_nr <- 
    ase_tlr |>
    group_split(donor_id, replic, var_id) |>
    future_map(quietly_run_model_norandom)

write_rds(tlr_res_nr, "./results_glm/tlr_nr.rds")

dn2_res_nr <- 
    ase_dn2 |>
    group_split(donor_id, replic, var_id) |>
    future_map(quietly_run_model_norandom)

write_rds(dn2_res_nr, "./results_glm/dn2_nr.rds")


##########

# Process results
read_results <- function(f) {
    read_rds(f) |>
    map("result") |>
    map_dfr("dat")
}


bcr_res_df <- read_results("results_glm/bcr.rds")
tlr_res_df <- read_results("results_glm/tlr.rds")
dn2_res_df <- read_results("results_glm/dn2.rds")
res_random_df <- bind_rows(bcr_res_df, tlr_res_df, dn2_res_df)

bcr_res_nr_df <- read_results("results_glm/bcr_nr.rds")
tlr_res_nr_df <- read_results("results_glm/tlr_nr.rds")
dn2_res_nr_df <- read_results("results_glm/dn2_nr.rds")
res_norand_df <- bind_rows(bcr_res_nr_df, tlr_res_nr_df, dn2_res_nr_df)

# Plot
library(patchwork)
library(extrafont)

stim_colors <- 
    c("Day 0" = "#898e9f",
      "BCR" = "#003967",
      "TLR7" = "#637b31",
      "DN2" = "#a82203")

qq_random_df <- 
    res_random_df |>
    group_split(donor_id, replic, stim) |>
    map_df(~arrange(., p) |> 
	   mutate(observed = -log10(p), 
		  expected = -log10(ppoints(n())))) |>
    unite("sample_id", c(donor_id, replic), sep = "_")

qq_norand_df <- 
    res_norand_df |>
    filter(!is.na(p)) |>
    group_split(donor_id, replic, stim) |>
    map_df(~arrange(., p) |> 
	   mutate(observed = -log10(p), 
		  expected = -log10(ppoints(n())))) |>
    unite("sample_id", c(donor_id, replic), sep = "_")


qq_plot <- 
    ggplot(qq_random_df) +
    geom_point(aes(x = expected, y = observed, color = stim),
	       size = 1, alpha = .5) +
    geom_abline(intercept = 0, slope = 1) +
    scale_color_manual(values = stim_colors) +
    facet_grid(sample_id~stim) +
    theme_minimal() +
    theme(legend.position = "none",
	  panel.grid.minor = element_blank(),
	  strip.text.y = element_text(angle = 0),
	  plot.background = element_rect(color = "white", fill = "white")) +
    labs(x = expression(paste("Expected -log"[10], plain(P))),
	 y = expression(paste("Observed -log"[10], plain(P))))

ggsave("./plots/qq_glm.png", qq_plot, height = 12, width = 4)

qq_nr_plot <- 
    ggplot(qq_norand_df) +
    geom_point(aes(x = expected, y = observed, color = stim),
	       size = 1, alpha = .5) +
    geom_abline(intercept = 0, slope = 1) +
    scale_color_manual(values = stim_colors) +
    facet_grid(sample_id~stim, scale = "free_y") +
    theme_minimal() +
    theme(legend.position = "none",
	  panel.grid.minor = element_blank(),
	  strip.text.y = element_text(angle = 0),
	  plot.background = element_rect(color = "white", fill = "white")) +
    labs(x = expression(paste("Expected -log"[10], plain(P))),
	 y = expression(paste("Observed -log"[10], plain(P))))

ggsave("./plots/qq_nr_glm.png", qq_nr_plot, height = 12, width = 4)


# Correlation between replicates
compute_cor <- function(s1, s2) {

    data_1 <- 
	filter(qq_norand_df, sample_id == s1) |>
	select(1:3, observed)
    
    data_2 <- 
	filter(qq_norand_df, sample_id == s2) |>
	select(1:3, observed)

   inner_join(data_1, data_2, join_by(stim, var_id)) |>
   group_by(stim) |>
   summarise(r = cor(observed.x, observed.y, method = "pearson"))
}

compute_pi1 <- function(s1, s2) {
    
    data_1 <- 
	filter(qq_norand_df, sample_id == s1) |>
	select(1:3, p) |>
	filter(qvalue::qvalue(p)$qvalue <= .2)
    
    data_2 <- 
	filter(qq_norand_df, sample_id == s2) |>
	select(1:3, p) |>
	semi_join(data_1, join_by(stim, var_id))

    data_2 |>
	group_by(stim) |>
	summarise(pi1 = 1L - qvalue::qvalue(p)$pi0)
}


plot_corr <- function(stim_i) {

    lim <- round(range(cor_df$r, na.rm = TRUE), 1)

    ggplot(cor_df |> filter(stim == stim_i), 
	   aes(x = sample_1, y = sample_2)) +
    geom_tile(aes(fill = r), alpha = .8) +
    scale_fill_gradient(low = "white", high = stim_colors[stim_i]) +
    geom_text(aes(label = round(r, 2)),
	      size = 4, fontface = "bold", family = "Arial") +
    scale_x_discrete(labels = function(x) sub("^(\\d+)_([ABC])$", "\\1 \\2", x))  +
    scale_y_discrete(labels = function(x) sub("^(\\d+)_([ABC])$", "\\1 \\2", x))  +
    theme_minimal() +
    theme(
	  axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1, family = "Arial"),
	  axis.text.y = element_text(size = 14, family = "Arial"),
	  plot.title = element_text(size = 16, hjust = .5, family = "Arial", face = "bold"),
	  panel.grid = element_blank(),
	  plot.background = element_rect(color = "white", fill = "white")) +
    labs(x = NULL, y = NULL, title = stim_i) +
    guides(fill = "none")
}

plot_pi1 <- function(stim_i) {

    lim <- round(range(pi_plot_df$pi1, na.rm = TRUE), 1)

    ggplot(pi_plot_df |> filter(stim == stim_i), 
	   aes(x = sample_1, y = sample_2)) +
    geom_tile(aes(fill = pi1), alpha = .8) +
    scale_fill_gradient(low = "white", high = stim_colors[stim_i]) +
    geom_text(aes(label = round(pi1, 2)),
	      size = 4, fontface = "bold", family = "Arial") +
    scale_x_discrete(labels = function(x) sub("^(\\d+)_([ABC])$", "\\1 \\2", x))  +
    scale_y_discrete(labels = function(x) sub("^(\\d+)_([ABC])$", "\\1 \\2", x))  +
    theme_minimal() +
    theme(
	  axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1, family = "Arial"),
	  axis.text.y = element_text(size = 14, family = "Arial"),
	  plot.title = element_text(size = 16, hjust = .5, family = "Arial", face = "bold"),
	  panel.grid = element_blank(),
	  plot.background = element_rect(color = "white", fill = "white")) +
    labs(x = NULL, y = NULL, title = stim_i) +
    guides(fill = "none")
}


pairwise_sample_combs <- 
    qq_norand_df |>
    distinct(sample_1 = sample_id) |>
    expand_grid(sample_2 = sample_1)

cor_df <-
    pairwise_sample_combs |>
    mutate(sample_1 = pmin(sample_1, sample_2),
	   sample_2 = pmax(sample_1, sample_2)) |>
    distinct(sample_1, sample_2) |>
    mutate(data = map2(sample_1, sample_2, compute_cor)) |>
    unnest(cols = data) |>
    arrange(sample_1, sample_2) |>
    mutate(sample_1 = fct_inorder(sample_1),
	   sample_2 = fct_rev(sample_2)) |>
    filter(sample_1 != sample_2)

corr_plot <- 
    plot_corr("BCR") + plot_corr("TLR7") + plot_corr("DN2") + plot_layout(nrow = 1)

ggsave("./plots/reps_corr_glm.png", corr_plot, width = 22, height = 7)

pi1_df <- 
    pairwise_sample_combs |>
    filter(sample_1 != sample_2) |>
    mutate(data = map2(sample_1, sample_2, compute_pi1)) |>
    unnest(cols = data)

pi_plot_df <- 
    pairwise_sample_combs |>
    filter(sample_1 == sample_2) |>
    mutate(pi1 = NA) |>
    bind_rows(pi1_df) |>
    arrange(sample_1, sample_2)

pi1_plot <- 
    plot_pi1("BCR") + plot_pi1("TLR7") + plot_pi1("DN2") + plot_layout(nrow = 1)

ggsave("./plots/reps_corr_glm.png", corr_plot, width = 22, height = 7)

# Plot betas between technical replicates
get_rep_betas <- function(d, s, r1, r2) {
    
    data_r1 <- 
	res_norand_df |>
	filter(donor_id == d, replic == r1, stim == s,
	       qvalue::qvalue(p)$qvalue < 0.1) |>
	select(donor_id, stim, var_id, beta_time)

    data_r2 <-
	res_norand_df |>
	filter(donor_id == d, replic == r2, stim == s) |>
	inner_join(data_r1, join_by(donor_id, stim, var_id)) |>
	select(var_id, beta_time.x, beta_time.y, p)

    data_r2
}

replic_data <- 
    res_norand_df |>
    distinct(donor_id = donor_id, replic, stim) |>
    group_by(donor_id, stim) |>
    filter(n_distinct(replic) > 1) |>
    ungroup()

replic_spec <-     
    left_join(replic_data, replic_data, 
	  join_by(donor_id, stim), 
	  relationship = "many-to-many") |>
    filter(replic.x != replic.y) |>
    select(donor_id, stim, replic.x, replic.y) |>
    arrange(donor_id, stim, replic.x, replic.y) |>
    mutate(data = pmap(list(donor_id, stim, replic.x, replic.y), get_rep_betas)) |>
    unnest(data)

plot_data <- 
    replic_spec |>
    mutate(lab = sprintf("Donor #%s\n Rep %s vs. Rep %s", donor_id, replic.x, replic.y))

cor_data <- 
    plot_data |>
    group_by(lab) |>
    mutate(ymax = max(beta_time.y)) |>
    group_by(stim) |>
    mutate(xmin = min(beta_time.x)) |>
    group_by(lab, donor_id, stim, replic.x, replic.y, xmin, ymax) |>
    summarise(r = cor(x = beta_time.x, y = beta_time.y, method = "spearman")) |>
    ungroup() |>
    mutate(cor_lab = paste("\u03C1 =", round(r, 2)))

p <- 
    ggplot(plot_data, aes(x = beta_time.x, y = beta_time.y)) +
    geom_abline() +
    geom_point(aes(color = stim), alpha = .75) +
    scale_color_manual(values = stim_colors) +
    scale_x_continuous(breaks = scales::pretty_breaks(3)) +
    scale_y_continuous(breaks = scales::pretty_breaks(3)) +
    geom_text(data = cor_data, aes(x = xmin, y = ymax + .25, label = cor_lab),
	      size = 3.5, family = "Arial", hjust = "inward") +
    facet_grid(lab~stim, scales = "free") +
    theme_minimal() +
    theme(legend.position = "none",
	  axis.title = element_blank(),
	  strip.text.x = element_text(size = 12, family = "arial", face = "bold"),
	  strip.text.y = element_text(size = 12, family = "arial", angle = 0),
	  panel.grid.minor = element_blank(),
	  plot.background = element_rect(fill = "white", color = "white"))

ggsave("./plots/betas_replic.png", p, height = 8.5, width = 5)

pvalues_plot <-
    plot_data |>
    mutate(p = replace_na(p, 1)) |>
    ggplot(aes(x = p, y = after_stat(density))) +
    geom_histogram(aes(fill = stim)) +
    scale_fill_manual(values = stim_colors) +
    scale_x_continuous(breaks = c(0, 1)) +
    facet_grid(lab~stim, scales = "free") +
    theme_minimal() +
    theme(legend.position = "none",
	  axis.text = element_text(size = 11, family = "Arial"),
	  axis.title = element_text(size = 11, family = "Arial"),
	  strip.text.x = element_text(size = 12, family = "arial", face = "bold"),
	  strip.text.y = element_text(size = 12, family = "arial", angle = 0),
	  panel.grid.minor = element_blank(),
	  plot.background = element_rect(fill = "white", color = "white")) +
    labs(x = "P-value")

ggsave("./plots/pvals_replic.png", pvalues_plot, height = 8.5, width = 5)


# Plot cases of significant effects of time

ase_data

signif_vars <- res_norand_df |> 
    filter(qvalue::qvalue(p)$qvalue < 0.05) |>
    arrange(p) |>
    distinct(var_id) |>
    pull(var_id)

ase_plot_data <- res_norand_df |>
    filter(var_id %in% signif_vars) |>
    mutate(var_id = factor(var_id, levels = signif_vars)) |>
    arrange(var_id, donor_id, replic, stim) |>
    inner_join(ase_data, join_by(donor_id, replic, stim, var_id)) |>
    mutate(stim = factor(stim, levels = names(stim_colors))) |>
    unite("sample_id", c(donor_id, replic), sep = "_")

ase_plot_data_counts <- ase_data |>
    pivot_longer(ref_count:alt_count, names_to = "allele", values_to = "counts") |>
    mutate(allele = str_remove(allele, "_count"))

p <- 
    ggplot(ase_plot_data_counts |> filter(var_id == first(var_id)), 
	   aes(x = counts, y = allele)) +
    geom_col(aes(fill = stim)) +
    #geom_text(data = ase_plot_data |> filter(var_id == first(var_id))) +
    scale_fill_manual(values = stim_colors) +
    facet_grid(sample_id~stim, scale = "free") +
    theme_minimal() +
    theme(plot.background = element_rect(fill = "white", color = "white"))

ggsave("./plots/ase_time.png", p)


