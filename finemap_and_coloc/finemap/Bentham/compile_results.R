library(tidyverse)

windows_df <- 
    read_tsv("./data/windows.tsv", col_names = c("locus", "coord")) |>
    rowid_to_column("window_id") |>
    mutate(window_id = as.character(window_id))

pip_df <- 
    sprintf("./data/susie/pip/window%s.tsv", windows_df$window_id) |>
    setNames(windows_df$window_id) |>
    map_dfr(data.table::fread, .id = "window_id") |>
    as_tibble() |>
    left_join(windows_df, join_by(window_id)) |>
    select(window_id, locus, rsid, ref, alt, pip, cs, coverage)

summ_stats <-
    read_tsv("./data/summary_stats.tsv") |>
    select(locus, rsid, logp)

# For coloc
pip_df |>
    filter(!is.na(cs)) |>
    left_join(summ_stats, join_by(locus, rsid)) |>
    filter(logp > -log10(5e-5)) |>
    distinct(window_id, locus, cs) |>
    arrange(as.integer(window_id), cs) |>
    print(n = Inf)

pip_df |> 
    filter(!is.na(cs)) |>
    arrange(as.integer(window_id), cs, desc(pip)) |>
    filter(locus == "TNFSF4")
    
"./data/susie/diagnostics/window3.tsv" |> 
    read_tsv() |>
    filter(grepl("rs2205960|rs10912578", snp_id))

pip_df |>
    left_join(summ_stats, join_by(locus, rsid)) |>
    filter(rsid %in% c("rs4728141", "rs2004640", "rs78658945"))


