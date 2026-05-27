library(tidyverse)
library(glue)
library(RColorBrewer)
library(pheatmap)

stims <- read_lines("./stims.txt")

homer_df <- 
    glue("results/{stims}/knownResults.txt") |>
    setNames(stims) |>
    map_dfr(~read_tsv(.) |> 
	    janitor::clean_names() |>
	    {function(x) setNames(x, str_remove(names(x), "_of_\\d+$"))}(),
	    .id = "stim")

plot_mat <- 
    homer_df |>
    mutate_at(vars(starts_with("percent")), parse_number) |>
    mutate(fc = percent_of_target_sequences_with_motif/percent_of_background_sequences_with_motif) |>
    group_by(motif_name) |>
    filter(any(q_value_benjamini <= 0.01),
	   any(percent_of_target_sequences_with_motif > 5),
	   any(fc > 1)) |>
    ungroup() |>
    select(stim, motif_name, fc) |>
    pivot_wider(names_from = stim, values_from = fc) |>
    column_to_rownames("motif_name") |>
    data.matrix()

png("./homer.png", units = "in", height = 25, width = 8.5, res = 300)
pheatmap(plot_mat, 
	 fontsize = 9, 
	 angle_col = 0,
	 scale = "row",
	 cluster_cols = FALSE,
	 color = colorRampPalette(c("white", "firebrick3"))(100))
dev.off()



dn2_annot_header <- 
    read_tsv("./results/DN2/annotate.txt", col_names = FALSE, n_max = 1) |>
    pivot_longer(everything()) |>
    mutate(value = case_when(grepl("^PeakID", value) ~ "PeakID",
			     TRUE ~ value)) |>
    pull(value)

dn2_annot <- 
    data.table::fread("./results/DN2/annotate.txt", skip = 1) |>
    as_tibble() |>
    setNames(dn2_annot_header) |>
    janitor::clean_names()

dn2_annot |>
    filter(peak_id == "Interval_134768") |>
    select(-(5:21)) |>
    pivot_longer(-(1:4)) |>
    filter(value != "") |>
    select(-(peak_id:end)) |>
    extract(value, c("dist_center", "info"), "(\\d+)(\\(.+\\))", convert = TRUE) |>
    arrange(dist_center) |>
    filter(grepl("^irf4|tbx21", name)) |>
    print(n = Inf)

dn2_annot |>
    filter(peak_id == "Interval_134768") |>
    select(start:end) |>
    mutate(center = (start + end)/2,
	   dist_var = 160028076 - center)

