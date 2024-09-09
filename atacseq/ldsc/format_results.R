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
ldsc_files <- 
    list.files("results", 
	       pattern = "cell_type_results.txt",
	       full.names = TRUE)

names(ldsc_files) <- 
    ldsc_files |> 
    basename() |>
    str_remove("^ldscseg_") |>
    str_remove("\\.cell_type_results.txt$")

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
    select(gwas, set = Name, tau_hat = Coefficient, p = Coefficient_P_value, pfdr, tau_star) |>
    ungroup()

write_tsv(out, "./compiled_results.tsv")
















#res <- 
#    read_lines("./module_results.log", skip = 29, n_max = 8) |>
#    {function(x) tibble(dat = x)}() |>
#    separate(dat, c("stat", "dat"), sep = ": ")
#
#annots <-
#    str_split(res$dat[[1]], " ") |>
#    unlist()
#
#res_table <- res |>
#    slice(-1) |>
#    separate(dat, annots, sep = "\\s+", convert = TRUE) |>
#    pivot_longer(-stat, names_to = "annot") |>
#    pivot_wider(names_from = stat) |>
#    janitor::clean_names()
#
#res_table |>
#    mutate(annot = recode(annot, 
#			  "L2_1" = "turquoise",
#			  "L2_2" = "blue",
#			  "L2_3" = "brown",
#			  "L2_4" = "yellow",
#			  "L2_5" = "green",
#			  "L2_6" = "red",
#			  "L2_7" = "black")) |>
#    select(annot, coefficients, coefficient_se) |>
#    mutate(z = coefficients/coefficient_se,
#	   p = pnorm(-abs(z)),
#	   p_adj = p.adjust(p, method = "fdr")) |>
#    tail(7)
#
#res250 <- 
#    read_tsv("./modules_250_results.results") |>
#    mutate(padj = p.adjust(Enrichment_p, method = "fdr")) |>
#    tail(7) |>
#    select(Category, Enrichment:padj) |>
#    print(width = Inf)
#
#res500 <- 
#    read_tsv("./module_500_results.results") |>
#    mutate(padj = p.adjust(Enrichment_p, method = "fdr")) |>
#    tail(7) |> 
#    select(Category, Enrichment:padj) |>
#    print(width = Inf)
#
#res250 |> select(Category, padj)
#res500 |> select(Category, padj)
#
#read_tsv("./module_500_results.results") |>
#    mutate(padj = p.adjust(Enrichment_p, method = "fdr")) |>
#    filter(padj <= 0.05) |>
#    arrange(padj) |> 
#    pull(Category)
