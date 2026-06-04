library(tidyverse)
library(data.table)
library(susieR)

# -----------------------------------------------------------------------------
# 1. LOAD LD REFERENCE DATA
# -----------------------------------------------------------------------------
sample_size <- 6748 + 11516

# Load the variant IDs from the PLINK2 .vars file
ld_vars <- 
    fread("./data/chr1_TNFSF4_r2.unphased.vcor1.vars", header = FALSE) |>
    pull(V1)

# Load the LD matrix
ld_matrix <- 
    fread("./data/chr1_TNFSF4_r2.unphased.vcor1") |> 
    data.matrix()

rownames(ld_matrix) <- colnames(ld_matrix) <- ld_vars

# -----------------------------------------------------------------------------
# 2. LOAD & MERGE GWAS DATA
# -----------------------------------------------------------------------------
# Read the VCF file
gwas_vcf_raw <- 
    read_tsv("./data/langefeld_gwas_grch37.vcf.gz", comment = "##")

gwas_final <- 
    gwas_vcf_raw |>
    select(snp_id = ID, FORMAT, stats = stats) |>
    separate_rows(FORMAT, stats, sep = ":") |>
    pivot_wider(names_from = FORMAT, values_from = stats) |>
    mutate(
	   ES = as.numeric(ES),
	   SE = as.numeric(SE),
	   LP = as.numeric(LP),
	   z = ES / SE
    ) |>
    # Keep only variants that match the LD panel
    filter(snp_id %in% colnames(ld_matrix))

# -----------------------------------------------------------------------------
# 3. SUBSET LD MATRIX
# -----------------------------------------------------------------------------
shared_snps <- gwas_final$snp_id
ld_final <- ld_matrix[shared_snps, shared_snps]

# -----------------------------------------------------------------------------
# 4. RUN SUSIE
# -----------------------------------------------------------------------------
condz <- kriging_rss(gwas_final$z, ld_final, n = sample_size, r_tol = 1e-04)

condz_df <- 
    condz[[2]] |>
    as_tibble() |>
    bind_cols(select(gwas_final, snp_id)) |>
    mutate(flag = ifelse(logLR > 2 & abs(z) > 2, TRUE, FALSE))

write_tsv(condz_df, "./results/susie_diagnostics_TNFSF4.tsv")

# Filter out the 9 SNPs flagged by the Kriging diagnostic
flagged_snps <- condz_df |> filter(flag == TRUE) |> pull(snp_id)

gwas_cleansed <- gwas_final |> filter(!snp_id %in% flagged_snps)
ld_cleansed <- ld_final[gwas_cleansed$snp_id, gwas_cleansed$snp_id]

fit <- 
    susie_rss(
	      z = gwas_cleansed$z, 
	      R = ld_cleansed, 
	      n = sample_size, 
	      L = 5, 
	      coverage = 0.9, 
	      r_tol = 1e-05
)

# -----------------------------------------------------------------------------
# 5. FORMAT & EXPORT RESULTS
# -----------------------------------------------------------------------------
if (!is.null(fit$sets$coverage)) {
    coverage_df <- tibble(
        coverage = scales::percent(round(fit$sets$coverage, 3)),
        cs = paste0("L", fit$sets$cs_index)
    )

    cs_df <- 
        enframe(fit$sets$cs, "cs", "rowid") |>
        unnest(cols = rowid) |>
        left_join(coverage_df, by = "cs")

    pip_df <- 
        enframe(fit$pip, "snp_id", "pip") |>
        rowid_to_column() |>
        left_join(cs_df, by = "rowid")
} else {
    pip_df <- 
        enframe(fit$pip, "snp_id", "pip") |>
        mutate(cs = NA, coverage = NA)
}

write_tsv(pip_df, "./results/susie_pip_TNFSF4.tsv")
