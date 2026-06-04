# ==============================================================================
# Description:  Aggregates the individual HOMER known motif enrichment results 
#               across all stimulation conditions into a single, clean master 
#               table. This is essential for downstream visualization (e.g., 
#               heatmaps, dot plots) and generating supplementary tables.
# Input:        1. ./data/stims.txt (List of conditions)
#               2. ./results_vert/{stims}/knownResults.txt (HOMER outputs)
# Output:       ./results_vert/results.tsv (Combined & cleaned master table)
# ==============================================================================

library(tidyverse)
library(glue)

# ------------------------------------------------------------------------------
# 1. Data Import & Aggregation
# ------------------------------------------------------------------------------
# Read the list of stimulation conditions generated during the setup step
stims <- read_lines("./data/stims.txt")

# Dynamically construct file paths, read them, and bind them into a single dataframe
results <-
    glue("./results/{stims}/knownResults.txt") |>
    setNames(paste0(stims, "c")) |>
    map_dfr(~read_tsv(.) |>
            janitor::clean_names() |>
            {function(x) setNames(x, str_remove(names(x), "_of_\\d+$"))}(),
            .id = "stim") |>
    mutate_at(vars(starts_with("percent_of_")), parse_number) |>
    mutate(stim = recode(stim, "IL4c" = "IL-4c"),
	   stim = factor(stim, levels = c("IL-4c", "TLR7c", "BCRc", "DN2c")),
           fc = percent_of_target_sequences_with_motif/percent_of_background_sequences_with_motif) |>
    extract(motif_name,
	    c("tf_name", "tf_family", "dataset", "db"), 
	    "(.+?)(?:\\((.+)\\))?/([^/]+)/(.+)", 
	    remove = FALSE) |>
    mutate(log10p = log_p_value/log(10)) |>
    select(stim, motif_name, tf_name, tf_family, dataset, consensus,
           pct_target =  percent_of_target_sequences_with_motif,
           pct_bg = percent_of_background_sequences_with_motif,
           fc, log10p, q_value = q_value_benjamini)

# ------------------------------------------------------------------------------
# 2. Export Table
# ------------------------------------------------------------------------------
write_tsv(results, "./results/results.tsv")

results |>
    filter(q_value < 0.01) |>
    write_tsv("./results/results_fdr01.tsv")

results |>
    filter(q_value < 0.01) |>
    distinct(motif_name) |>
    write_tsv("./results/motifs_fdr01.tsv")

