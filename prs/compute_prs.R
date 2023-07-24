library(tidyverse)
library(readxl)
library(ggbeeswarm)
library(furrr)

sle_df <-
    "./mgb_data/sle_data.tsv" |>
    read_tsv(col_types = c(.default = "c"))

summ_stats <- 
    "./common_data/Khunsriraksakul_summstats_hg38.tsv" |>
    read_tsv()

vcf_df <- 
    "./common_data/allchr.mgb.pt_5e-4.vcf.gz" |>
    data.table::fread() |>
    as_tibble() |>
    select(-(REF:FORMAT)) |>
    rename(chr = `#CHROM`) |>
    left_join(select(summ_stats, chr, pos, rsid, beta), 
	      join_by(chr, POS == pos)) |>
    select(chr, POS, ID, rsid, beta, everything())

plan(multisession, workers = availableCores())

dose_df <- 
    vcf_df |>
    group_split(ID) |>
    future_map_dfr(~pivot_longer(., -(chr:beta), names_to = "sample_id", values_to = "genot") |>
		   extract("sample_id", "subject_id", "[^-]+-(\\d+)") |>
		   separate(genot, c("a1", "a2"), sep = "\\|", convert = TRUE) |>
		   mutate(dose = map2_int(a1, a2, sum)) |>
		   select(chr, POS, ID, rsid, beta, subject_id, dose))

prs_df <- 
    dose_df |>
    mutate(upt_dose = ifelse(beta < 0, dose, 2L - dose),
	   upt_beta = ifelse(beta < 0, beta * -1L, beta),
	   rs = upt_dose * upt_beta) |>
    group_by(subject_id) |>
    summarise(prs = sum(rs)) |>
    ungroup()

plot_df <- inner_join(sle_df, prs_df, join_by(subject_id))
write_tsv(plot_df, "./common_data/prs_p+t_5e-4.tsv")

p <- ggplot(plot_df, aes(x = group, y = prs, fill = group)) +
    geom_violin(alpha = .5) +
    scale_fill_manual(values = c("Control" = "midnightblue", "SLE" = "tomato4")) +
    theme_bw() +
    theme(legend.position = "none",
	  plot.title = element_text(size = 8)) +
    labs(x = NULL,
	 y = "PRS",
	 title = expression(paste("PRS computed on 828 variants from P+T 5x", 10^{-4}, "method")))


ggsave("./plots/prs.png", p, width = 4, height = 4)
