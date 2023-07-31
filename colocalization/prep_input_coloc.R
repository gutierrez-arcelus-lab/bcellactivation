library(tidyverse)
library(rvest)

# functions
download_tbi <- function(ftp) {
    
    ftptbi <- paste0(ftp, ".tbi")
    idx_file <- file.path(getwd(), basename(ftptbi))

    if (! file.exists(idx_file) ) download.file(ftptbi, dest = idx_file)
}


# eQTL catalogue tabix data
tabix_paths <- 
    "https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv" |>
    read_tsv()

imported_tabix_paths <- 
    "https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths_imported.tsv" |>
    read_tsv()

paths_df <- bind_rows(select(tabix_paths, -sample_size), imported_tabix_paths)

paths_df |>
    select(study, tissue_label, condition_label, qtl_group, quant_method, ftp_path) |>
    write_tsv("./data/coloc_input/eqtl_catalogue_paths.tsv")

paths_df |>
    filter(study != "GTEx_V8") |>
    pull(ftp_path) |>
    write_lines("./data/coloc_input/ftp_urls.txt")

# save column names
read_tsv(paths_df$ftp_path[1], n_max = 1) |>
    colnames() |>
    write_lines("./eqtl_column_names.txt")

# download index files
# instead of downloading again, I used soft links
#walk(paths_df$ftp_path, download_tbi)

# GWAS summ stats
# in both the GWAS catalog and OpenGWAS data there are variants with same rsid and different stats 
# From OpenGWAS
gwas_stats <- "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Bentham/ebi-a-GCST003156.vcf.gz" |>
    read_tsv(comment = "##") |>
    mutate(`#CHROM` = paste0("chr", `#CHROM`)) |>
    add_count(`#CHROM`, ID) |>
    filter(n == 1) |>
    select(-n)

# Harmonized data from GWAS Catalog
#gwas_stats <- "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Bentham/bentham_GRCh38.tsv.gz" |>
#    read_tsv() |>
#    drop_na(hm_odds_ratio) |>
#    select(chr = hm_chrom, rsid = hm_rsid, pos = hm_pos, 
#	   eff_allele = hm_effect_allele, other_allele = hm_other_allele,
#	   odds_ratio, beta, se = standard_error, p = p_value) |>
#    mutate(chr = paste0("chr", chr))

# liftOver
chain <- "/reference_databases/ReferenceGenome/liftover_chain/hg19/hg19ToHg38.over.chain.gz"
bed19 <- "./data/gwas/Bentham_OpenGWAS_hg19.bed"  
bed38 <- "./data/gwas/Bentham_OpenGWAS_hg38.bed" 
fail <- "./data/gwas/Bentham_OpenGWAS_failToLift.txt"

gwas_stats_hg19 <- gwas_stats |>
    select(chr = `#CHROM`, start = POS, rsid = ID) |>
    mutate(end = start, 
	   start = start - 1) |>
    select(chr, start, end, rsid)

write_tsv(gwas_stats_hg19, bed19, col_names = FALSE)

"liftOver %s %s %s %s" |>
    sprintf(bed19, chain, bed38, fail) |>
    system()

gwas_stats_38 <- bed38 |>
    read_tsv(col_names = c("chr", "start", "end", "rsid")) |>
    select(chr, pos = end, rsid) |>
    left_join(gwas_stats, by = c("chr" = "#CHROM", "rsid" = "ID")) |>
    select(chr, pos, rsid, ref = REF, alt = ALT, stats_info = FORMAT, stats = `EBI-a-GCST003156`)

# Top variants
bentham_top <- "https://www.nature.com/articles/ng.3434/tables/1" |>
    read_html() |>
    html_node("table") |>
    html_table(header = TRUE, fill = TRUE) |>
    select(rsid = 1, chr = 2, pos = 3, locus = 4, p = 5, or = 6) |>
    slice(-1) |>
    as_tibble() |>
    mutate(rsid = gsub("\\s|[a-zA-Z,]+$", "", rsid),
           pos = parse_number(pos),
           p = gsub("\\s", "", p),
	   p = sub("^([0-9.]+)[^0-9.]10[^0-9.](\\d+)$", "\\1e-\\2", x = p),
	   p = parse_number(p),
           or = parse_number(or),
	   chr = paste0("chr", chr))

# Save a description of the region
top_regions <- bentham_top |>
    filter(p < 1e-5) |>
    select(chr, rsid, locus) |>
    inner_join(gwas_stats_38) |>
    select(chr, pos, locus) |>
    mutate(chr = factor(chr, levels = c(unique(chr), numeric = TRUE))) |>
    arrange(chr, pos) |>
    rowid_to_column("region")

write_tsv(top_regions, "./data/bentham_top_hits.tsv")

tier2_regions <- bentham_top |>
    filter(between(p, 1e-5, 5e-3)) |>
    select(chr, rsid, locus) |>
    inner_join(gwas_stats_38) |>
    select(chr, pos, locus) |>
    mutate(chr = factor(chr, levels = c(unique(chr), numeric = TRUE))) |>
    arrange(chr, pos) |>
    rowid_to_column("region") |>
    mutate(region = region + max(top_regions$region))

write_tsv(tier2_regions, "./data/bentham_tier2_hits.tsv")

bentham_top_pos <- bentham_top |>
    filter(p < 1e-5) |>
    select(chr, rsid) |>
    inner_join(gwas_stats_38) |>
    select(chr, pos) |>
    mutate(chr = factor(chr, levels = c(unique(chr), numeric = TRUE))) |>
    arrange(chr, pos) |>
    rowid_to_column("region")

bentham_tier2_pos <- bentham_top |>
    filter(between(p, 1e-5, 5e-3)) |>
    select(chr, rsid) |>
    inner_join(gwas_stats_38) |>
    select(chr, pos) |>
    mutate(chr = factor(chr, levels = c(unique(chr), numeric = TRUE))) |>
    arrange(chr, pos) |>
    rowid_to_column("region") |>
    mutate(region = region + max(top_regions$region))

coords <- bentham_top_pos |>
    mutate(start = pos - 5e5, end = pos + 5e5,
	   chr = sub("chr", "", chr),
	   coord = sprintf("%s:%d-%d", chr, start, end))

coords_tier2 <- bentham_tier2_pos |>
    mutate(start = pos - 5e5, end = pos + 5e5,
	   chr = sub("chr", "", chr),
	   coord = sprintf("%s:%d-%d", chr, start, end))

# Save region for eQTL catalogue query
coords |>
    select(region, coord) |>
    write_tsv("./data/coloc_input/regions_bentham.tsv")

coords_tier2 |>
    select(region, coord) |>
    write_tsv("./data/coloc_input/regions_bentham_tier2.tsv")



# Write minimum GWAS dataset for coloc
gwas_min_df <- bentham_top_pos |>
    mutate(data = map2(chr, pos, ~filter(gwas_stats_38, chr == .x, between(pos, .y - 5e5, .y + 5e5)))) |>
    select(region, data) |>
    unnest(cols = data) |>
    separate(stats, c("beta", "se", "log10p", "dummy"), sep = ":", convert = TRUE) |>
    mutate(gwas_varbeta = se^2) |>
    select(region, chr, pos, rsid, gwas_beta = beta, gwas_varbeta)

write_tsv(gwas_min_df, "./data/coloc_input/gwas_data_bentham.tsv")

gwas_min_tier2_df <- bentham_tier2_pos |>
    mutate(data = map2(chr, pos, ~filter(gwas_stats_38, chr == .x, between(pos, .y - 5e5, .y + 5e5)))) |>
    select(region, data) |>
    unnest(cols = data) |>
    separate(stats, c("beta", "se", "log10p", "dummy"), sep = ":", convert = TRUE) |>
    mutate(gwas_varbeta = se^2) |>
    select(region, chr, pos, rsid, gwas_beta = beta, gwas_varbeta)

write_tsv(gwas_min_tier2_df, "./data/coloc_input/gwas_data_bentham_tier2.tsv")


# Select genes 
annotations <- 
    file.path("/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens", 
              "gencode.v30.primary_assembly.annotation.gtf.gz") |>
    read_tsv(comment = "#", col_types = "c-cii-c-c",
             col_names = c("chr", "feature", "start", "end", "strand", "info"))

eqtl_regions <- bentham_top_pos |>
    mutate(left = pos - 2.5e5, right = pos + 2.5e5)

# increase distance
eqtl_regions_tier2 <- bentham_tier2_pos |>
    mutate(left = pos - 4e5, right = pos + 4e5)

bed <- annotations |>
    filter(feature == "gene") |>
    mutate(tss = case_when(strand == "+" ~ start, 
			   strand == "-" ~ end, 
			   TRUE ~ NA_integer_)) |>
    inner_join(eqtl_regions, join_by(chr, between(tss, left, right))) |>
    mutate(gene_id = str_extract(info, "(?<=gene_id\\s\")[^\"]+"),
	   gene_name = str_extract(info, "(?<=gene_name\\s\")[^\"]+"),
	   gene_id = sub("^(ENSG\\d+)\\.\\d+((_PAR_Y)?)$", "\\1\\2", gene_id)) |>
    select(region, chr, start, end, gene_id, gene_name)

bed_tier2 <- annotations |>
    filter(feature == "gene") |>
    mutate(tss = case_when(strand == "+" ~ start, 
			   strand == "-" ~ end, 
			   TRUE ~ NA_integer_)) |>
    inner_join(eqtl_regions_tier2, join_by(chr, between(tss, left, right))) |>
    mutate(gene_id = str_extract(info, "(?<=gene_id\\s\")[^\"]+"),
	   gene_name = str_extract(info, "(?<=gene_name\\s\")[^\"]+"),
	   gene_id = sub("^(ENSG\\d+)\\.\\d+((_PAR_Y)?)$", "\\1\\2", gene_id)) |>
    select(region, chr, start, end, gene_id, gene_name)

gene_ids <- select(bed, region, gene_id, gene_name)
gene_ids_tier2 <- select(bed_tier2, region, gene_id, gene_name)

# Save gene IDs for eQTL catalogue filtering
write_tsv(gene_ids, "./data/coloc_input/gene_annots_bentham.tsv")
write_tsv(gene_ids_tier2, "./data/coloc_input/gene_annots_bentham_tier2.tsv")


# Langefeld
langefeld_tier1 <- "../sle_variants/paper_data/langefeld_tableS2.xlxs" |>
    readxl::read_excel(1, skip = 2) |>
    select(snp_id = 1, 
           chr = 2, 
           pos = 3, 
           gene_region = 4, 
           region_rank = 5, 
           ref_allele = 6, 
           p = `P-value`,
           or = `OR (95% CI)`) |>
    filter(!is.na(pos)) |>
    mutate(snp_id = gsub("\\s|[a-zA-Z,]+$", "", snp_id),  
           p = gsub("\\s", "", p),
	   p = sub("^([0-9.]+)[^0-9.]10[^0-9.](\\d+)$", "\\1e-\\2", p),
	   p = parse_number(p),
           or = parse_number(or)) |>
    select(snp_id, chr, pos, gene_region, region_rank, ref_allele, p, or)

make_lange_min_df <- function(dat, chrom) {

    dat |>
    filter(grepl("^rs", rsid)) |>
    extract(rsid, "rsid", "(rs\\d+)") |>
    select(rsid, position, beta = `frequentist_add_beta_1:add/sle=1`, 
	   se = frequentist_add_se_1) |>
    drop_na() |>
    add_count(rsid) |>
    filter(n == 1) |>
    mutate(varbeta = se^2) |>
    select(rsid, position, beta, varbeta) |>
    inner_join(filter(langefeld_top_pos, chr == chrom), join_by(between(position, left, right))) |>
    select(region, chr, pos = position, rsid, gwas_beta = beta, gwas_varbeta = varbeta) |> 
    left_join(snp_ids_map, by = c("chr", "rsid" = "rsid")) |>
    left_join(snp_ids_map, by = c("chr", "rsid" = "rsid_new")) |>
    filter(!is.na(pos_hg38.x) | !is.na(pos_hg38.y)) |>
    mutate(final_rsid = case_when(!is.na(rsid_new) & rsid == rsid_new ~ rsid,
				  !is.na(rsid_new) & rsid != rsid_new ~ rsid_new,
				  is.na(rsid_new) ~ rsid.y,
				  TRUE ~ NA_character_),
	   pos_hg38 = case_when(!is.na(pos_hg38.x) ~ pos_hg38.y,
				!is.na(pos_hg38.y) ~ pos_hg38.x,
				TRUE ~ NA_integer_)) |>
    select(region, chr, pos = pos_hg38, rsid = final_rsid, gwas_beta, gwas_varbeta)
}

langefeld_top_pos <- langefeld_tier1 |>
    filter(grepl("IL10|IKZF1", gene_region)) |>
    mutate(chr = sub("^(\\d+)[pq]\\d+$", "\\1", chr)) |>
    select(region = region_rank, chr, pos) |>
    mutate(left = pos - 5e5, right = pos + 5e5)

snp_ids_map <- "./data/snpout.txt.result.txt" |>
    read_tsv(col_names = c("chr", "rsid", "pos_hg38", "rsid_new"), col_types = "ccic")


langefeld_1 <- "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Langefeld/ea.imputed.chr1.out" |>
    read_delim(delim = " ", comment = "#")

langefeld_7 <- "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Langefeld/ea.imputed.chr7.out" |>
    read_delim(delim = " ", comment = "#")

min_df_lange_1 <- make_lange_min_df(langefeld_1, chrom = 1)
min_df_lange_7 <- make_lange_min_df(langefeld_7, chrom = 7)

min_df_lange <- bind_rows(min_df_lange_1, min_df_lange_7)

write_tsv(min_df_lange, "./data/coloc_input/gwas_data_langefeld.tsv")

langefeld_top_pos_hg38 <- langefeld_tier1 |>
    filter(grepl("IL10|IKZF1", gene_region)) |>
    mutate(chr = sub("^(\\d+)[pq]\\d+$", "\\1", chr)) |>
    select(region = region_rank, chr, rsid = snp_id) |>
    left_join(snp_ids_map, by = c("chr", "rsid")) |>
    select(region, chr, rsid, pos = pos_hg38) |> 
    mutate(left = pos - 5e5, 
	   right = pos + 5e5)

langefeld_top_pos_hg38 |>
    mutate(coord = paste0(chr, ":", left, "-", right)) |>
    select(region, coord) |>
    write_tsv("./data/coloc_input/regions_langefeld.tsv")

# Save gene IDs for eQTL catalogue filtering
bed_langefeld <- annotations |>
    filter(feature == "gene") |>
    mutate(tss = case_when(strand == "+" ~ start, 
			   strand == "-" ~ end, 
			   TRUE ~ NA_integer_)) |>
    inner_join(langefeld_top_pos_hg38 |> mutate(chr = paste0("chr", chr)), 
	       join_by(chr, between(tss, left, right))) |>
    mutate(gene_id = str_extract(info, "(?<=gene_id\\s\")[^\"]+"),
	   gene_name = str_extract(info, "(?<=gene_name\\s\")[^\"]+"),
	   gene_id = sub("^(ENSG\\d+)\\.\\d+((_PAR_Y)?)$", "\\1\\2", gene_id)) |>
    select(chr, start, end, gene_id, gene_name)

bed_langefeld |>
    select(gene_id, gene_name) |>
    write_tsv("./data/coloc_input/gene_annots_langefeld.tsv")

