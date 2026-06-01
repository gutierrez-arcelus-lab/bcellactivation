# ==============================================================================
# Description:  Aggregates the raw LDSC outputs across all GWAS traits into a 
#               single master table. It extracts the annotation standard deviation 
#               and total heritability from the S-LDSC models, merges them with 
#               the LDSC-SEG coefficients, applies FDR correction, and calculates 
#               tau-star ($\tau^*$) for standardized cross-trait comparisons.
# Input:        1. 1000G baseline .M_5_50 files (SNP counts)
#               2. ./data/traits.txt
#               3. ./results/sldsc_*.results & *.log (S-LDSC outputs)
#               4. ./results/ldscseg_*.cell_type_results.txt (LDSC-SEG outputs)
# Output:       ./compiled_results.tsv (Master heritability enrichment table)
# ==============================================================================

library(tidyverse)
library(glue)

# ------------------------------------------------------------------------------
# 1. Compute M (Total SNPs) from Baseline Reference
# ------------------------------------------------------------------------------
# M is the total number of SNPs in the baseline reference panel used for LDSC.
m_h2 <- 
    paste0("./data/ldsc_data/v4/1000G_Phase3_baseline_v1.2_ldscores/baseline.", 1:22, ".l2.M_5_50") |>
    setNames(1:22) |>
    map_dfr(~read_tsv(., col_names = FALSE) |> select(X1), .id = "chrom") |>
    summarise(m = sum(X1)) |>
    pull(m)

# ------------------------------------------------------------------------------
# 2. Extract Data from S-LDSC Models
# ------------------------------------------------------------------------------
# Load the GWAS trait metadata
traits <- 
    read_tsv("./data/traits.txt", col_names = c("directory", "trait", "gwas", "ref"))

# A. Calculate Annotation Standard Deviation (SD)
# Extracts the proportion of SNPs (prop_snps) falling into our ATAC peaks from 
# the S-LDSC results files. For a binary annotation, SD = sqrt(p(1-p)).
annot_df <- 
    glue("./results/sldsc_{traits$gwas}.results") |>
    setNames(traits$gwas) |>
    map_dfr(~read_tsv(.) |>
	    janitor::clean_names() |>
	    filter(grepl("^L2", category)) |>
	    select(annot = category, prop_snps = prop_sn_ps),
	    .id = "trait") |>
    mutate(annot = str_remove(annot, "L2"),
	   annot = str_remove(annot, "_\\d$"),
	   sdev = sqrt(prop_snps * (1 - prop_snps)))

# B. Extract Total Heritability
# Parses the raw S-LDSC log files to extract the "Total Observed scale h2" value.
h2_df <- 
    glue("./results/sldsc_{traits$gwas}.log") |>
    setNames(traits$gwas) |>
    map_dfr(~read_lines(.) |>
	    {function(x) keep(x, grepl("^Total Observed scale h2", x))}() |>
	    str_extract("h2: ([0-9.]+)", group = TRUE) |>
	    as.numeric() |>
	    {function(x) tibble(h2 = x)}(),
	    .id = "trait")
 
# ------------------------------------------------------------------------------
# 3. Extract LDSC-SEG Results & Compute tau*
# ------------------------------------------------------------------------------

# Parse the main LDSC-SEG coefficient and p-value tables
ldsc_df <- 
    glue("results/ldscseg_{traits$gwas}.cell_type_results.txt") |>
    read_tsv(id = "trait") |>
    mutate(trait = str_remove(trait, "results/ldscseg_"),
	   trait = str_remove(trait, ".cell_type_results.txt"))

# Merge all extracted metrics, apply FDR multiple testing correction, and 
# calculate the standardized effect size (tau*). 
# tau* represents the proportionate change in per-SNP heritability associated 
# with a 1 standard deviation increase in the value of the annotation, allowing 
# direct comparison across traits with different total heritabilities.
out <- 
    ldsc_df |>
    left_join(annot_df, join_by(trait, Name == annot)) |>
    left_join(h2_df, join_by(trait)) |>
    mutate(pfdr = p.adjust(Coefficient_P_value, method = "fdr")) |>
    rowwise() |>
    mutate(tau_star = ( (m_h2 * sdev)/h2 ) * Coefficient) |>
    select(trait, set = Name, tau_hat = Coefficient, p = Coefficient_P_value, pfdr, tau_star) |>
    ungroup()

# Export final compiled results
write_tsv(out, "./compiled_results.tsv")
