library(tidyverse)
library(readxl)
library(glue)

parse_pvalues <- function(p) {

    p |>
	{function(x) gsub("\\s", "", x)}() |>
	{function(x) sub("^([0-9.]+)[^0-9.]10[^0-9.](\\d+)$", "\\1e-\\2", x)}() |>
        parse_number()
}

read_langefeld <- function(tab_path, sheet_number) {
    read_excel(tab_path, sheet_number, skip = 2) |>
    select(snp_id = 1, 
           chr = 2, 
           pos = 3, 
           gene_region = 4, 
           region_rank = 5, 
           ref_allele = 6, 
           p = `P-value`,
           or = `OR (95% CI)`,
           p_stepwise = `Regional Stepwise  P-value`) |> 
    filter(!is.na(pos)) |>
    mutate(snp_id = gsub("\\s|[a-zA-Z,]+$", "", snp_id),  
           p = parse_pvalues(p),
           p_stepwise = parse_pvalues(p_stepwise),
           or = parse_number(or)) |>
    select(snp_id, chr, pos, gene_region, region_rank, ref_allele, p, p_stepwise, or)
}

read_summ_stats <- function(chrom) { 
    glue("/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Langefeld/ea.imputed.chr{chrom}.out") |>
    data.table::fread(skip = "alternate_ids") |>
    as_tibble() |>
    filter(all_maf >= 0.005) |>
    select(rsid, 
	   pos = position, 
	   alleleA, 
	   alleleB,
	   beta = `frequentist_add_beta_1:add/sle=1`, 
	   se = frequentist_add_se_1, 
	   p = frequentist_add_wald_pvalue_1) |>
    drop_na(beta, se) |>
    mutate(chr = chrom) |>
    select(chr, everything())
}

# Create directories to save results
if (! file.exists("data")) dir.create("data")
if (! file.exists("data/ld")) dir.create("data/ld")
if (! file.exists("data/susie")) dir.create("data/susie")
if (! file.exists("data/susie/diagnostics")) dir.create("data/susie/diagnostics")
if (! file.exists("data/susie/log")) dir.create("data/susie/log")
if (! file.exists("data/susie/pip")) dir.create("data/susie/pip")
if (! file.exists("data/susie/lbf")) dir.create("data/susie/lbf")

tier1_ea <-
    "../../../sle_variants/paper_data/langefeld_tableS2.xlxs" |>
    read_langefeld(1) |>
    group_by(gene_region) |>
    slice(which.min(p)) |>
    ungroup() |>
    mutate(chr = sub("^(\\d+)[qp]\\d+$", "\\1", chr)) |>
    mutate(chr = factor(chr, levels = c(1:22, "X", "Y"))) |>
    arrange(chr, pos)

write_tsv(tier1_ea, "./data/sentinels.tsv")

tier1_regions <- 
    tier1_ea |>
    mutate(start = pos - 5e5, end = pos + 5e5) |>
    mutate(region = glue("{chr}:{start}-{end}")) |>
    select(gene_region, chr, start, end, region)

tier1_regions |>
    select(gene_region, region) |>
    write_tsv("./data/windows.tsv", col_names = FALSE)

gwas_chromosomes <- as.character(unique(tier1_regions$chr))

gwas_stats_regions <- 
    gwas_chromosomes |>
    map_dfr(~read_summ_stats(.) |>
	    inner_join(tier1_regions, join_by(chr, between(pos, start, end))) |>
	    select(chr, gene_region, rsid, pos, alleleA, alleleB, beta, se, p)) |>
    mutate(rsid = str_extract(rsid, "rs\\d+"))

# Match variants to 1000G data
kgp_data <- 
    "/reference_databases/1000G_VCF/phase3/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz" |>
    read_tsv(comment = "##") |> 
    select(-QUAL, -FILTER, -INFO) |>
    inner_join(tier1_regions, join_by(`#CHROM` == chr, between(POS, start, end))) |> 

    select(chr = `#CHROM`, pos = POS, rsid = ID, ref = REF, alt = ALT) |> 
    separate_rows(alt, sep = ",") |> 
    distinct()

variant_stats_final <- 
    gwas_stats_regions |>
    distinct(chr, rsid, pos, alleleA, alleleB) |>
    inner_join(kgp_data, join_by(chr, pos, alleleA == ref, alleleB == alt)) |>
    select(chr, rsid = rsid.y, pos, alleleA, alleleB)

gwas_stats_final <-
    gwas_stats_regions |>
    select(-rsid) |>
    left_join(variant_stats_final, join_by(chr, pos, alleleA, alleleB)) |>
    select(gene_region, chr, rsid, pos, alleleA, alleleB, beta, se, p)

# Save GWAS data
write_tsv(gwas_stats_final, "./data/summary_stats.tsv")


# Save KGP IDs for LD reference panel
dat <- 
    "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index" |>
    read_tsv(comment = "##")

europeans <- c("CEU", "TSI", "GBR", "IBS", "FIN")

dat |>
    distinct(SAMPLE_NAME, POPULATION) |>
    filter(POPULATION %in% europeans) |>
    arrange(SAMPLE_NAME) |>
    pull(SAMPLE_NAME) |>
    write_lines("./data/ref_panel.txt")
