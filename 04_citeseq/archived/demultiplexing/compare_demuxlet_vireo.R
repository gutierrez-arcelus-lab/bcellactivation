library(tidyverse)

donors <- read_lines("./genotypes/data/mgb_sample_ids.txt")

demuxlet_df <- 
    "./demuxlet/results/demuxlet_calls.tsv" |>
    read_tsv()

vireo_df <-
    "./vireo/results/vireo_calls.tsv" |>
    read_tsv()

demux_df <-
    full_join(demuxlet_df, vireo_df, 
	      join_by(barcode, batch), suffix = c("_demuxlet", "_vireo"))

demux_tally <- 
    demux_df |>
    count(donor_id_demuxlet, donor_id_vireo, sort = TRUE) |>
    mutate(donor_id_demuxlet = factor(donor_id_demuxlet, levels = c(donors, "DBL", "AMB")),
	   donor_id_vireo = factor(donor_id_vireo, levels = c(donors, "doublet", "unassigned"))) |>
    complete(donor_id_demuxlet, donor_id_vireo, fill = list(n = 0))

p <- 
    ggplot(demux_tally, aes(x = donor_id_demuxlet, y = donor_id_vireo)) +
    geom_tile(aes(fill = n)) +
    geom_text(aes(label = n), size = 8, size.unit = "pt", color = "black") +
    scale_fill_gradient(low = "white", high = "tomato") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    labs(x = "Demuxlet", y = "vireo") +
    guides(fill = "none")

ggsave("testp.png", p, width = 4.25, height = 4)

# Compare with HTODemux
htodemux_df <- 
    "../htodemux_calls.tsv" |>
    read_tsv() |>
    select(barcode, batch = orig.ident, hto_call = hto_global)

demux_df |>
    filter(status_vireo == "unassigned") |>
    select(barcode, batch, vireo_call = donor_id_vireo, demuxlet_call = donor_id_demuxlet) |>
    left_join(htodemux_df, join_by(barcode, batch)) |>
    count(vireo_call, demuxlet_call, hto_call, sort = TRUE) |>
    print(n = Inf)
