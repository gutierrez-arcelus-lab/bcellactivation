library(tidyverse)
library(glue)

# Compute M SNPs from Reference
m_h2 <- 
    paste0("/lab-share/IM-Gutierrez-e2/Public/tools/ldsc_data/v4/1000G_Phase3_baseline_v1.2_ldscores/baseline.", 1:22, ".l2.M_5_50") |>
    setNames(1:22) |>
    map_dfr(~read_tsv(., col_names = FALSE) |> select(X1), .id = "chrom") |>
    summarise(m = sum(X1)) |>
    pull(m)

# Traits
traits <- 
    read_tsv("./data/traits.txt", col_names = c("directory", "trait", "gwas", "ref"))

# Annotation SD
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

# h2
h2_df <- 
    glue("./results/sldsc_{traits$gwas}.log") |>
    setNames(traits$gwas) |>
    map_dfr(~read_lines(.) |>
	    {function(x) keep(x, grepl("^Total Observed scale h2", x))}() |>
	    str_extract("h2: ([0-9.]+)", group = TRUE) |>
	    as.numeric() |>
	    {function(x) tibble(h2 = x)}(),
	    .id = "trait")
  
# Compile LDSC-SEG results
ldsc_df <- 
    glue("results/ldscseg_{traits$gwas}.cell_type_results.txt") |>
    read_tsv(id = "trait") |>
    mutate(trait = str_remove(trait, "results/ldscseg_"),
	   trait = str_remove(trait, ".cell_type_results.txt"))

out <- 
    ldsc_df |>
    left_join(annot_df, join_by(trait, Name == annot)) |>
    left_join(h2_df, join_by(trait)) |>
    mutate(pfdr = p.adjust(Coefficient_P_value, method = "fdr")) |>
    rowwise() |>
    mutate(tau_star = ( (m_h2 * sdev)/h2 ) * Coefficient) |>
    select(trait, set = Name, tau_hat = Coefficient, p = Coefficient_P_value, pfdr, tau_star) |>
    ungroup()

write_tsv(out, "./compiled_results.tsv")
