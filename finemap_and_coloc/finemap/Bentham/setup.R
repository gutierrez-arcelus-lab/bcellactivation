library(tidyverse)
library(rvest)

if (! file.exists("data")) dir.create("data")
if (! file.exists("data/ld")) dir.create("data/ld")
if (! file.exists("data/susie")) dir.create("data/susie")
if (! file.exists("data/susie/diagnostics")) dir.create("data/susie/diagnostics")
if (! file.exists("data/susie/log")) dir.create("data/susie/log")
if (! file.exists("data/susie/pip")) dir.create("data/susie/pip")
if (! file.exists("data/susie/lbf")) dir.create("data/susie/lbf")


# Summary statistics
gwas_stats <- 
    "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Bentham/opengwas/ebi-a-GCST003156.vcf.gz" |>
    read_tsv(comment = "##") |>
    select(chrom = `#CHROM`, pos = POS, rsid = ID, ref = REF, alt = ALT, FORMAT, stats = starts_with("EBI"))

# Bentham et al lead variants
bentham_top <- 
    "https://www.nature.com/articles/ng.3434/tables/1" |>
    read_html() |>
    html_node("table") |>
    html_table(header = TRUE, fill = TRUE) |>
    select(rsid = 1, chrom = 2, pos = 3, locus = 4, p = 5, or = 6) |>
    slice(-1) |>
    as_tibble() |>
    mutate(rsid = gsub("\\s|[a-zA-Z,]+$", "", rsid),
           pos = parse_number(pos),
           p = gsub("\\s", "", p),
	   p = sub("^([0-9.]+)[^0-9.]10[^0-9.](\\d+)$", "\\1e-\\2", p),
	   p = parse_number(p),
           or = parse_number(or)) |> 
    filter(p <= 5e-5) |>
    filter(locus != "MHC class IIId")

write_tsv(bentham_top, "./data/sentinels.tsv")

# Use ±500kb windows around top variants
bentham_top_windows <- 
    bentham_top |>
    mutate(start = pos - 5e5L, end = pos + 5e5L,
	   locus = gsub(" ", "", locus)) |>
    select(locus, chrom, start, end)

gwas_stats_windows <- 
    inner_join(x = gwas_stats, y = bentham_top_windows, 
	       join_by(chrom, between(pos, start, end))) |>
    select(locus, chrom, pos, rsid, ref, alt, FORMAT, stats) |>
    separate_rows(FORMAT:stats, sep = ":") |>
    pivot_wider(names_from = FORMAT, values_from = stats) |>
    select(locus:alt, beta = ES, se = SE, logp = LP)

write_tsv(gwas_stats_windows, "./data/summary_stats.tsv")


# Save windows
windows_coords <- 
    gwas_stats_windows |>
    mutate(locus = fct_inorder(locus)) |>
    group_by(locus) |>
    summarise(coord = sprintf("%s:%d-%d", unique(chrom), min(pos), max(pos))) |>
    ungroup()

write_tsv(windows_coords, "./data/windows.tsv", col_names = FALSE)


# Reference panel for LD
kgp_df <- 
    "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index" |>
    read_tsv(comment = "##")

kgp_df |>
    distinct(SAMPLE_NAME, POPULATION) |>
    filter(POPULATION %in% c("CEU", "TSI", "GBR", "FIN", "IBS")) |>
    arrange(SAMPLE_NAME) |>
    pull(SAMPLE_NAME) |>
    write_lines("./data/ref_panel.txt")
