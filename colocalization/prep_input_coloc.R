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

# download index files
#walk(paths_df$ftp_path, download_tbi)


# GWAS summ stats
# From OpenGWAS
gwas_stats <- "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Bentham/ebi-a-GCST003156.vcf.gz" |>
    read_tsv(comment = "##") |>
    mutate(`#CHROM` = paste0("chr", `#CHROM`))

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

#gwas_stats <- "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Bentham/bentham_GRCh38.tsv.gz" |>
#    read_tsv() |>
#    drop_na(hm_odds_ratio) |>
#    select(chr = hm_chrom, rsid = hm_rsid, pos = hm_pos, 
#	   eff_allele = hm_effect_allele, other_allele = hm_other_allele,
#	   odds_ratio, beta, se = standard_error, p = p_value) |>
#    mutate(chr = paste0("chr", chr))

# Top variants
bentham_top <- "https://www.nature.com/articles/ng.3434/tables/1" |>
    read_html() |>
    html_node("table") |>
    html_table(header = TRUE, fill = TRUE) |>
    select(rsid = 1, chr = 2, pos = 3, locus = 4, p = 11, or = 12) |>
    slice(-1) |>
    as_tibble() |>
    mutate(rsid = gsub("\\s|[a-zA-Z,]+$", "", rsid),
           pos = parse_number(pos),
           p = gsub("\\s", "", p) |>
	    sub("^([0-9.]+)[^0-9.]10[^0-9.](\\d+)$", "\\1e-\\2", x = _) |>
	    parse_number(),
           or = parse_number(or),
	   chr = paste0("chr", chr))

bentham_top_pos <- bentham_top |>
    select(chr, rsid) |>
    inner_join(gwas_stats_38) |>
    select(chr, pos)

coords <- bentham_top_pos |>
    mutate(start = pos - 5e5, end = pos + 5e5,
	   chr = sub("chr", "", chr),
	   region = sprintf("%s:%d-%d", chr, start, end))

# Save region for eQTL catalogue query
write_lines(coords$region, "./data/coloc_input/region.txt")


# Write minimum GWAS dataset for coloc
gwas_min_df <- bentham_top_pos |>
    rowid_to_column("region") |>
    mutate(data = map2(chr, pos, ~filter(gwas_stats_38, chr == .x, between(pos, .y - 5e5, .y + 5e5)))) |>
    select(region, data) |>
    unnest(cols = data) |>
    separate(stats, c("beta", "se", "log10p", "dummy"), sep = ":", convert = TRUE) |>
    mutate(gwas_varbeta = se^2) |>
    select(region, chr, pos, rsid, gwas_beta = beta, gwas_varbeta)

write_tsv(gwas_min_df, "./data/coloc_input/gwas_data.tsv")


# Select genes 
annotations <- 
    file.path("/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens", 
              "gencode.v30.primary_assembly.annotation.gtf.gz") |>
    read_tsv(comment = "#", col_types = "c-cii-c-c",
             col_names = c("chr", "feature", "start", "end", "strand", "info"))

eqtl_regions <- bentham_top_pos |>
    mutate(left = pos - 2.5e5, right = pos + 2.5e5)

bed <- annotations |>
    filter(feature == "gene") |>
    mutate(tss = case_when(strand == "+" ~ start, 
			   strand == "-" ~ end, 
			   TRUE ~ NA_integer_)) |>
    inner_join(eqtl_regions, join_by(chr, between(tss, left, right))) |>
    mutate(gene_id = str_extract(info, "(?<=gene_id\\s\")[^\"]+"),
	   gene_name = str_extract(info, "(?<=gene_name\\s\")[^\"]+"),
	   gene_id = sub("^(ENSG\\d+)\\.\\d+((_PAR_Y)?)$", "\\1\\2", gene_id)) |>
    select(chr, start, end, gene_id, gene_name)

gene_ids <- select(bed, gene_id, gene_name)

# Save gene IDs for eQTL catalogue filtering
write_tsv(gene_ids, "./data/coloc_input/gene_annots.tsv")


