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

if (!file.exists("data")) dir.create("data")

tier1_ea <-
    "../sle_variants/paper_data/langefeld_tableS2.xlxs" |>
    read_langefeld(1) |>
    group_by(gene_region) |>
    slice_min(p) |>
    filter(n() == 1 | (n() > 1 & is.na(p_stepwise))) |>
    ungroup() |>
    mutate(chr = sub("^(\\d+)[qp]\\d+$", "\\1", chr)) |>
    mutate(chr = factor(chr, levels = c(1:22, "X", "Y"))) |>
    arrange(chr, pos)

write_tsv(tier1_ea, "./data/langefeld_sentinels.tsv")


tier1_regions <- 
    tier1_ea |>
    mutate(start = pos - 2.5e5, end = pos + 2.5e5) |>
    mutate(region = glue("{chr}:{start}-{end}")) |>
    select(gene_region, chr, start, end, region)

tier1_regions |>
    select(region, gene_region) |>
    write_tsv("./data/langefeld_regions.tsv", col_names = FALSE)

gwas_chromosomes <- as.character(unique(tier1_regions$chr))

gwas_stats_regions <- 
    gwas_chromosomes |>
    map_dfr(~read_summ_stats(.) |>
	    inner_join(tier1_regions, join_by(chr, between(pos, start, end))) |>
	    select(chr, gene_region, rsid, pos, alleleA, alleleB, beta, se, p)) |>
    mutate(rsid = str_extract(rsid, "rs\\d+"))


# Update and get missing RsIDs by matching position and alleles to dbSNP
dbsnp_vcf <- 
    "/lab-share/IM-Gutierrez-e2/Public/References/dbSNP/GCF_000001405.25.gz"

dbsnp_meta <- 
    "/lab-share/IM-Gutierrez-e2/Public/References/dbSNP/GCF_000001405.25_GRCh37.p13_assembly_report.txt" |>
    data.table::fread(skip = "# Sequence-Name") |>
    as_tibble() |>
    select(chr = `UCSC-style-name`, chr_accn = `RefSeq-Accn`) |>
    filter(chr %in% glue("chr{gwas_chromosomes}"))

dbsnp_cmds <- 
    tier1_regions |>
    select(gene_region, chr, start, end) |>
    mutate(chr = paste0("chr", chr)) |>
    left_join(dbsnp_meta) |>
    mutate(region = glue("{chr_accn}:{start}-{end}")) |>
    select(gene = gene_region, region) |>
    mutate(awk = "| awk '{ print $2,$3,$4,$5 }'",
	   cmd = glue("tabix {dbsnp_vcf} {region} {awk} > data/dbsnp_{gene}.vcf")) |>
    pull(cmd)

walk(dbsnp_cmds, system)

dbsnp_data <- 
    tier1_regions |>
    mutate(dbsnp_file = glue("data/dbsnp_{gene_region}.vcf")) |>
    select(gene_region, dbsnp_file) |>
    deframe() |>
    map_dfr(~read_delim(., delim = " ", col_names = c("pos", "rsid", "ref", "alt")) |>
	    separate_rows(alt, sep = ","), .id = "gene_region")
    
gwas_stats_regions_2 <- 
    gwas_stats_regions |>
    left_join(dbsnp_data, join_by(gene_region, pos, alleleA == ref, alleleB == alt)) |>
    filter(!is.na(rsid.x), !is.na(rsid.y), rsid.x != rsid.y)
    mutate(rsid = ifelse(!is.na(rsid.y), rsid.y, rsid.x)) |>
    select(chr, gene_region, rsid, pos, alleleA, alleleB, beta, se, p) |>
    filter(!is.na(rsid))

# Save GWAS data
write_tsv(gwas_stats_regions_2, "./data/langefeld_summ_stats.tsv")




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
    write_lines("./data/eur_kgp.txt")




















# Save coords in hg19 for LD estimation
paste0("2:", stat4_risk_var - 1e6, "-", stat4_risk_var + 1e6) |>
    write_lines("./data/stat4_coords_grch37.txt")

# liftOver Langefeld data
# LiftOver
bed19_file <- "./data/langefeld_hg19.bed"
bed38_file <- "./data/langefeld_hg38.bed"
fail_file <- "./data/langefeld.failTolift.txt"
chain_file <- "/reference_databases/ReferenceGenome/liftover_chain/hg19/hg19ToHg38.over.chain.gz"

bed19 <- 
    "./data/langefeld_stat4.tsv" |>
    read_tsv() |>
    mutate(chr = "chr2",
	   end = pos,
	   start = pos - 1L) |>
    select(chr, start, end, rsid)

write_tsv(bed19, bed19_file, col_names = FALSE)

command <- sprintf("liftOver %s %s %s %s", bed19_file, chain_file, bed38_file, fail_file)
system(command)

