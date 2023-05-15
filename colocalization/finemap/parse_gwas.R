library(tidyverse)
library(rvest)

# GWAS summ stats
gwas_stats <- "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Bentham/ebi-a-GCST003156.vcf.gz" |>
    read_tsv(comment = "##") |>
    mutate(`#CHROM` = paste0("chr", `#CHROM`)) |>
    add_count(`#CHROM`, ID) |>
    filter(n == 1) |>
    select(-n)

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
	   p = sub("^([0-9.]+)[^0-9.]10[^0-9.](\\d+)$", "\\1e-\\2", p),
	   p = parse_number(p),
           or = parse_number(or),
	   chr = paste0("chr", chr))

# liftOver
chain <- "/reference_databases/ReferenceGenome/liftover_chain/hg19/hg19ToHg38.over.chain.gz"
bed19 <- "./data/bentham_opengwas_hg19.bed"  
bed38 <- "./data/bentham_opengwas_hg38.bed" 
fail <- "./data/bentham_opengwas_failtolift.txt"

stats_hg19 <- gwas_stats |>
    select(chr = `#CHROM`, end = POS, varid = ID) |>
    mutate(start = end - 1L) |>
    select(chr, start, end, varid)

write_tsv(stats_hg19, bed19, col_names = FALSE)

"liftOver %s %s %s %s" |>
    sprintf(bed19, chain, bed38, fail) |>
    system()

stats_hg38 <- bed38 |>
    read_tsv(col_names = c("chr", "start", "end", "varid")) |>
    select(chr, bp_hg38 = end, varid) |>
    left_join(gwas_stats, join_by(chr == `#CHROM`, varid == ID)) |>
    select(chr, varid, bp = bp_hg38, ref = REF, alt = ALT, info = FORMAT, 
	   stats = `EBI-a-GCST003156`)

# Extract variants in windows around lead variants
lead_vars_windows <- bentham_top |>
    select(chr, varid = rsid, locus) |>
    inner_join(stats_hg38) |>
    mutate(window_start = bp - 2.5e5, 
	   window_end = bp + 2.5e5) |>
    select(chr, locus, window_start, window_end)

stats_windows <- 
    inner_join(stats_hg38, lead_vars_windows, 
	       join_by(chr, between(bp, window_start, window_end))) |>
    select(chr, locus, varid, bp, ref, alt, info, stats) |>
    separate(stats, c("beta", "se", "logp", "id"), sep = ":") |>
    select(-info, -id)

# Save regions
lead_vars_windows |>
    mutate(coord = sprintf("%s:%d-%d", chr, window_start, window_end)) |>
    pull(coord) |>
    write_lines("./data/regions_bentham.txt")




## Write minimum GWAS dataset for coloc
#gwas_min_df <- bentham_top_pos |>
#    mutate(data = map2(chr, pos, ~filter(gwas_stats_38, chr == .x, between(pos, .y - 5e5, .y + 5e5)))) |>
#    select(region, data) |>
#    unnest(cols = data) |>
#    separate(stats, c("beta", "se", "log10p", "dummy"), sep = ":", convert = TRUE) |>
#    mutate(gwas_varbeta = se^2) |>
#    select(region, chr, pos, rsid, gwas_beta = beta, gwas_varbeta)
#
#write_tsv(gwas_min_df, "./data/coloc_input/gwas_data_bentham.tsv")
#
## Select genes 
#annotations <- 
#    file.path("/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens", 
#              "gencode.v30.primary_assembly.annotation.gtf.gz") |>
#    read_tsv(comment = "#", col_types = "c-cii-c-c",
#             col_names = c("chr", "feature", "start", "end", "strand", "info"))
#
#eqtl_regions <- bentham_top_pos |>
#    mutate(left = pos - 2.5e5, right = pos + 2.5e5)
#
#bed <- annotations |>
#    filter(feature == "gene") |>
#    mutate(tss = case_when(strand == "+" ~ start, 
#			   strand == "-" ~ end, 
#			   TRUE ~ NA_integer_)) |>
#    inner_join(eqtl_regions, join_by(chr, between(tss, left, right))) |>
#    mutate(gene_id = str_extract(info, "(?<=gene_id\\s\")[^\"]+"),
#	   gene_name = str_extract(info, "(?<=gene_name\\s\")[^\"]+"),
#	   gene_id = sub("^(ENSG\\d+)\\.\\d+((_PAR_Y)?)$", "\\1\\2", gene_id)) |>
#    select(region, chr, start, end, gene_id, gene_name)
#
#gene_ids <- select(bed, region, gene_id, gene_name)
#
## Save gene IDs for eQTL catalogue filtering
#write_tsv(gene_ids, "./data/coloc_input/gene_annots_bentham.tsv")
#
