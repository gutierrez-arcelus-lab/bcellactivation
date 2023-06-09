library(tidyverse)
library(rvest)

# Bentham et al ###############################################################

# Harmonized GWAS summ stats from OpenGWAS
gwas_stats <- "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Bentham/ebi-a-GCST003156.vcf.gz" |>
    read_tsv(comment = "##") |>
    mutate(`#CHROM` = paste0("chr", `#CHROM`)) |>
    add_count(`#CHROM`, ID) |>
    filter(n == 1) |>
    select(-n)

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

# Use ±500kb windows around top variants instead of LD-based windows
bentham_top_windows <- bentham_top |>
    mutate(start = pos - 5e5L, stop = pos + 5e5L,
	   locus = gsub(" ", "", locus))

gwas_stats_windows <- 
    inner_join(x = gwas_stats, y = bentham_top_windows, 
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
###############################################################################

# Langefeld et al
langefeld_top <- read_tsv("../../sle_variants/paper_data/langefeld_top.tsv") |>
    mutate(chr = paste0("chr", chr))

langefeld_regions <- langefeld_top |>
    mutate(start = pos - 5e5L, end = pos + 5e5L) |>
    select(locus, chr, start, end)


langefeld_stats <- 
    "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Langefeld/ea.imputed.chr%d.out" |>
    sprintf(1:22) |>
    setNames(sprintf("chr%d", 1:22)) |>
    map_dfr(~read_delim(., comment = "#", delim = " ") |>
	    select(varid = rsid, bp = position, major = alleleA, minor = alleleB,
		   beta = `frequentist_add_beta_1:add/sle=1`, se = frequentist_add_se_1, 
		   p = frequentist_add_wald_pvalue_1) |>
	    drop_na(beta) |> 
	    filter(grepl("^rs", varid)) |>
	    extract(varid, "varid", "(rs\\d+)") |>
	    add_count(varid) |>
	    filter(n == 1) |>
	    select(-n), .id = "chr")

langefeld_stats_windows <- 
    inner_join(langefeld_stats, langefeld_regions, 
	       join_by(chr, between(bp, start, end))) |>
    select(chr, locus, varid, bp, major, minor, beta, se, p)

langefeld_bed19_file <- "./data/langefeld_1MbWindows_hg19.bed"  
langefeld_bed38_file <- "./data/langefeld_1MbWindows_hg38.bed" 
langefeld_fail <- "./data/langefeld_1MbWindows_failtolift.txt"

langefeld_bed19 <- langefeld_stats_windows |>
    mutate(start = bp - 1L) |>
    select(chr, start, end = bp, varid, locus)

write_tsv(langefeld_bed19, langefeld_bed19_file, col_names = FALSE)
    
"liftOver %s %s %s %s" |>
    sprintf(langefeld_bed19_file, chain, langefeld_bed38_file, langefeld_fail) |>
    system()
    
langefeld_stats_hg38 <- langefeld_bed38_file |>
    read_tsv(col_names = c("chr", "start", "end", "varid", "locus")) |>
    select(chr, bp = end, varid, locus) |>
    left_join(langefeld_stats_windows, join_by(chr, locus, varid)) |>
    select(chr, locus, varid, bp = bp.x, major, minor, beta, se, p) |>
    filter(chr %in% paste0("chr", c(1:22, "X")))

write_tsv(langefeld_stats_hg38, "./data/langefeld_1MbWindows_hg38_summstats.tsv")

# Save regions
langefeld_regions_df <- langefeld_stats_hg38 |>
    mutate(locus = fct_inorder(locus)) |>
    group_by(locus) |>
    summarise(coord = sprintf("%s:%d-%d", unique(chr), min(bp), max(bp))) |>
    ungroup()

langefeld_regions_df |>
    select(coord, locus) |>
    write_tsv("./data/regions_langefeld_1mb.tsv", col_names = FALSE)

# most major/minor alleles seem to match REF/ALT
for (i in langefeld_regions_df$locus) { 
    
    variants <- langefeld_stats_hg38 |>
	filter(locus == i) |>
        unite(ID, c(chr, bp, major, minor), sep = ":") |>
	pull(ID)

    out <- langefeld_regions_df |>
	filter(locus == i) |>
	pull(coord) |>
	{function(x) sprintf("./data/%s_langefeld.variants.txt", x)}()

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

