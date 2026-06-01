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
    # Append 'c' to the stimulation names (e.g., 'BCR' -> 'BCRc') to match the 
    # specific B cell state nomenclature used in the manuscript text.
    setNames(paste0(stims, "c")) |>
    # Iterate over the files, read the TSVs, and append a 'stim' identifier column
    map_dfr(~read_tsv(.) |>
	    # Standardize column names (removes spaces, special characters)
	    janitor::clean_names() |>
	    # HOMER appends dynamic numbers to some column names (e.g., "_of_450").
            # This regex strips those numeric suffixes so the columns align 
            # perfectly across all conditions when binding the rows together.
	    {function(x) setNames(x, str_remove(names(x), "_of_\\d+$"))}(),
	    .id = "stim")

# ------------------------------------------------------------------------------
# 2. Export Table
# ------------------------------------------------------------------------------
write_tsv(results, "./results/results.tsv")
