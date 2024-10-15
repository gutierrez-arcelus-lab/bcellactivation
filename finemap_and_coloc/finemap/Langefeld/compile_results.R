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

write_tsv(pip_df, "./data/compiled_susie_pips.tsv")
