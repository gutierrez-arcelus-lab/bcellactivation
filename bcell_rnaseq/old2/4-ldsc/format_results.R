library(tidyverse)
library(glue)

# Function to normalize coefficient (tau*)
parse_log <- function(trait) {

    logfile <-
	glue("./results/sldsc_{trait}.log") |>
	read_lines()

    h2 <- 
	keep(logfile, grepl("^Total Observed scale h2", logfile)) |>
	str_extract("h2: ([0-9.]+)", group = TRUE) |>
	as.numeric()

    M <- 
	keep(logfile, grepl("\\d+ SNPs remain", logfile)) |>
	last() |>
	str_extract("\\((\\d+) SNPs remain", group = TRUE) |>
	as.integer()

    tibble("h2" = h2, "M" = M)
}


# Compile LDSC-SEG results
traits <- 
    read_tsv("./data/traits.tsv", col_names = c("directory", "trait", "gwas", "ref"))

ldsc_files <-
    glue("results/ldscseg_{traits$gwas}.cell_type_results.txt") |>
    setNames(traits$gwas)

ldsc_df <- 
    map_dfr(ldsc_files, read_tsv, .id = "gwas")

annots <- 
    expand_grid(set = unique(ldsc_df$Name), chr = 1:22) |>
    mutate(f = glue("./data/ldscores/{set}.{chr}.annot"),
	   annot = map(f, ~read_tsv(.) |> pull(ANNOT)))

annot_sd <-
    annots |>
    select(set, annot) |>
    unnest(annot) |>
    group_by(set) |>
    summarise(sdev = sd(annot)) |>
    ungroup()

ldsc_estimates <- 
    ldsc_df |>
    distinct(gwas) |>
    mutate(data = map(gwas, parse_log)) |>
    unnest("data")

out <- 
    ldsc_df |>
    left_join(annot_sd, join_by(Name == set)) |>
    left_join(ldsc_estimates, join_by(gwas)) |>
    mutate(pfdr = p.adjust(Coefficient_P_value, method = "fdr")) |>
    rowwise() |>
    mutate(tau_star = ( (M * sdev)/h2 ) * Coefficient) |>
    select(gwas, set = Name, tau_hat = Coefficient, se = Coefficient_std_error, p = Coefficient_P_value, pfdr, tau_star) |>
    ungroup()

write_tsv(out, "./compiled_results.tsv")

