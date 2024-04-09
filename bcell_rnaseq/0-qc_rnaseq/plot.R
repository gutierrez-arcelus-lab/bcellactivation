library(tidyverse)
library(glue)
library(tidytext)

meta <- read_tsv("./metadata.tsv")

fastqs <- 
    meta |>
    select(donor_id, replic_id, stim) |>
    unite("sample_id", c(donor_id, replic_id, stim), sep = "_", remove = FALSE) |>
    mutate(fastq = glue("/temp_work/ch229163/fastq/highinput_merged/{donor_id}_{replic_id}_{stim}_R1.fastq.gz"))

total_qc_reads <- 
    "/temp_work/ch229163/fastq/highinput_merged/multiqc_data/multiqc_general_stats.txt" |>
    read_tsv() |> 
    select(sample_id = Sample, passed_reads = ends_with("total_sequences")) |>
    drop_na(passed_reads) |>
    separate(sample_id, c("donor_id", "replic_id", "stim", "r", "dummy1", "dummy2"), sep = "_") |>
    filter(r == "R1") |>
    select(donor_id, replic_id, stim, passed_reads)

total_reads <- 
    fastqs |>
    select(sample_id, fastq) |>
    deframe() |>
    map_dfr(~. |> 
	    paste0("_trimming_report.txt") |> 
	    read_lines() |> 
	    {function(x) grep("Total reads processed:", x, value = TRUE)}() |>
	    parse_number(), 
	    .id = "sample_id") |>
    pivot_longer(everything(), names_to = "sample_id", values_to = "reads") |>
    separate(sample_id, c("donor_id", "replic_id", "stim"), sep = "_")

reads_df <- 
    left_join(total_reads, total_qc_reads, join_by(donor_id, replic_id, stim)) |>
    select(donor_id, replic_id, stim, total = reads, passed = passed_reads) |>
    mutate(total = total - passed) |>
    pivot_longer(total:passed, names_to = "read_type", values_to = "n") |>
    unite("sample_id", c(donor_id, replic_id), sep = "_") |>
    mutate(stim = factor(stim, levels = c("unstday0", "BCR", "TLR7", "DN2")),
	   read_type = factor(read_type, levels = c("total", "passed")))

fill_colors <- 
    c("unstday0" = "grey",
      "BCR" = "royalblue",
      "TLR7" = "forestgreen",
      "DN2" = "tomato3")

p <- 
    ggplot(reads_df, 
       aes(x = n, 
	   y = reorder_within(sample_id, by = n, within = stim),
	   fill = stim,
	   alpha = read_type)) +
    geom_col() +
    geom_vline(xintercept = 40e6, linetype = 2, linewidth = .5) +
    scale_x_continuous(breaks = seq(0, 80e6, 20e6),
		       labels = function(x) x/1e6L) +
    scale_y_reordered() +
    scale_alpha_manual(values = c("total" = 1, "passed" = .7)) +
    scale_fill_manual(values = fill_colors) +
    facet_wrap(~stim, ncol = 2, scales = "free_y") + 
    theme_minimal() +
    theme(panel.grid = element_blank(),
	  plot.background = element_rect(color = "white", fill = "white")) +
    labs(x = "Million reads", y = NULL)

ggsave("./plots/total_reads.png", p)
