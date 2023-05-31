library(tidyverse)
library(rvest)

# Harmonized GWAS summ stats from OpenGWAS
gwas_stats <- "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Bentham/ebi-a-GCST003156.vcf.gz" |>
    read_tsv(comment = "##") |>
    mutate(`#CHROM` = paste0("chr", `#CHROM`)) |>
    add_count(`#CHROM`, ID) |>
    filter(n == 1) |>
    select(-n)

# LD blocks as used by Kundu et al, Nature Genetics, 2022
#ld_bed <- 
#    "https://bitbucket.org/nygcresearch/ldetect-data/raw/ac125e47bf7ff3e90be31f278a7b6a61daaba0dc/EUR/fourier_ls-all.bed" |>
#    read_tsv()
#
#write_tsv(ld_bed, "./data/fourier_ls-all.bed")

# At this step, variants at chrX are removed because there is no LD info for chrX.
# Anyway chrX would be removed later on, because there are not much data in the summ stats;
# Need to check if OpenGWAS harmonization is removing the data, 
# or if that already is the case in the original data.
#
#ld_bed <- read_tsv("./data/fourier_ls-all.bed") |>
#    rowid_to_column("ld_window")
#
# Bentham et al lead variants
#bentham_top <- "https://www.nature.com/articles/ng.3434/tables/1" |>
#    read_html() |>
#    html_node("table") |>
#    html_table(header = TRUE, fill = TRUE) |>
#    select(rsid = 1, chr = 2, pos = 3, locus = 4, p = 5, or = 6) |>
#    slice(-1) |>
#    as_tibble() |>
#    mutate(rsid = gsub("\\s|[a-zA-Z,]+$", "", rsid),
#           pos = parse_number(pos),
#           p = gsub("\\s", "", p),
#	   p = sub("^([0-9.]+)[^0-9.]10[^0-9.](\\d+)$", "\\1e-\\2", p),
#	   p = parse_number(p),
#           or = parse_number(or),
#	   chr = paste0("chr", chr)) |>
#    inner_join(ld_bed, join_by(chr, between(pos, start, stop)))
#    
#write_tsv(bentham_top, "./data/bentham_leadvars.tsv")

bentham_top <- read_tsv("./data/bentham_leadvars.tsv")

# liftOver
#chain <- "/reference_databases/ReferenceGenome/liftover_chain/hg19/hg19ToHg38.over.chain.gz"
#bed19 <- "./data/bentham_opengwas_ldwindows_hg19.bed"  
#bed38 <- "./data/bentham_opengwas_ldwindows_hg38.bed" 
#fail <- "./data/bentham_opengwas_failtolift.txt"

#stats_hg19 <- gwas_stats |>
#    select(chr = `#CHROM`, pos = POS, varid = ID) |>
#    inner_join(select(bentham_top, chr, ld_window, start, stop),
#	       join_by(chr, between(pos, start, stop))) |>
#    mutate(start = pos - 1L) |>
#    select(chr, start, end = pos, varid, ld_window)
#
#write_tsv(stats_hg19, bed19, col_names = FALSE)
#
#"liftOver %s %s %s %s" |>
#    sprintf(bed19, chain, bed38, fail) |>
#    system()
#
#stats_hg38 <- bed38 |>
#    read_tsv(col_names = c("chr", "start", "end", "varid", "ld_window")) |>
#    select(chr, bp = end, varid, ld_window) |>
#    left_join(gwas_stats, join_by(chr == `#CHROM`, varid == ID)) |>
#    left_join(select(bentham_top, locus, ld_window), join_by(ld_window)) |>
#    select(chr, ld_window, locus, varid, bp, ref = REF, alt = ALT, info = FORMAT, 
#	   stats = `EBI-a-GCST003156`) |>
#    separate_rows(c(info, stats), sep = ":", convert = TRUE) |>
#    pivot_wider(names_from = info, values_from = stats) |>
#    select(chr, ld_window, locus, varid, bp, ref, alt, beta = ES, se = SE, logp = LP) |>
#    mutate_at(vars(beta, se, logp), as.numeric)
#
#write_tsv(stats_hg38, "./data/bentham_opengwas_ldwindows_hg38_summstats.tsv")
#stats_hg38 <- read_tsv("./data/bentham_opengwas_ldwindows_hg38_summstats.tsv")
#
## Save regions
#regions_df <- stats_hg38 |>
#    group_by(ld_window, locus) |>
#    summarise(coord = sprintf("%s:%d-%d", unique(chr), min(bp), max(bp))) |>
#    ungroup()
#
#regions_df |>
#    select(coord, locus) |>
#    write_tsv("./data/regions_bentham.tsv", col_names = FALSE)
#
## Save variants
#for (i in regions_df$ld_window) { 
#    
#    variants <- stats_hg38 |>
#	filter(ld_window == i) |>
#        unite(ID, c(chr, bp, ref, alt), sep = ":") |>
#	pull(ID)
#
#    out <- regions_df |>
#	filter(ld_window == i) |>
#	pull(coord) |>
#	{function(x) sprintf("./data/%s.variants.txt", x)}()
#
#    write_lines(variants, out)
#}
#
#
#
# Use ±500kb windows around top variants instead of LD-based windows
bentham_top <- bentham_top |>
    mutate(start = pos - 5e5L, stop = pos + 5e5L,
	   locus = gsub(" ", "", locus))

gwas_stats_windows <- 
    inner_join(x = gwas_stats, y = bentham_top, 
	       join_by(`#CHROM` == chr, between(POS, start, stop))) |>
    select(locus, chr = `#CHROM`, POS, ID, REF, ALT, FORMAT, stats = starts_with("EBI"))

chain <- "/reference_databases/ReferenceGenome/liftover_chain/hg19/hg19ToHg38.over.chain.gz"
bed19_file <- "./data/bentham_opengwas_1MbWindows_hg19.bed"  
bed38_file <- "./data/bentham_opengwas_1MbWindows_hg38.bed" 
fail <- "./data/bentham_opengwas_1MbWindows_failtolift.txt"

bed19 <- gwas_stats_windows |>
    mutate(start = POS - 1L) |>
    select(chr, start, end = POS, ID, locus)

write_tsv(bed19, bed19_file, col_names = FALSE)
    
"liftOver %s %s %s %s" |>
    sprintf(bed19_file, chain, bed38_file, fail) |>
    system()

gwas_stats_windows_hg38 <- bed38_file |>
    read_tsv(col_names = c("chr", "start", "end", "varid", "locus")) |>
    select(chr, bp = end, varid, locus) |>
    left_join(gwas_stats_windows, join_by(chr, locus, varid == ID)) |>
    select(chr, locus, varid, bp, ref = REF, alt = ALT, info = FORMAT, stats) |>
    separate_rows(c(info, stats), sep = ":", convert = TRUE) |>
    pivot_wider(names_from = info, values_from = stats) |>
    select(chr, locus, varid, bp, ref, alt, beta = ES, se = SE, logp = LP) |>
    mutate_at(vars(beta, se, logp), as.numeric)

write_tsv(gwas_stats_windows_hg38, "./data/bentham_opengwas_1MbWindows_hg38_summstats.tsv")

# Save regions
regions_df <- gwas_stats_windows_hg38 |>
    mutate(locus = fct_inorder(locus)) |>
    group_by(locus) |>
    summarise(coord = sprintf("%s:%d-%d", unique(chr), min(bp), max(bp))) |>
    ungroup()

regions_df |>
    select(coord, locus) |>
    write_tsv("./data/regions_bentham_1mb.tsv", col_names = FALSE)

# Save variants
for (i in regions_df$locus) { 
    
    variants <- gwas_stats_windows_hg38 |>
	filter(locus == i) |>
        unite(ID, c(chr, bp, ref, alt), sep = ":") |>
	pull(ID)

    out <- regions_df |>
	filter(locus == i) |>
	pull(coord) |>
	{function(x) sprintf("./data/%s.variants.txt", x)}()

    write_lines(variants, out)
}


# Save IDs of Europeans in 1000G
dat <- "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index" |>
    read_tsv(comment = "##")

europeans <- c("CEU", "TSI", "GBR", "FIN", "IBS")

dat |>
    distinct(SAMPLE_NAME, POPULATION) |>
    filter(POPULATION %in% europeans) |>
    arrange(SAMPLE_NAME) |>
    pull(SAMPLE_NAME) |>
    write_lines("./data/europeans_samples.txt")

