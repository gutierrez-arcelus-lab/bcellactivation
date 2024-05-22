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

if (!file.exists("data")) dir.create("data")


tier1_ea <-
    "../sle_variants/paper_data/langefeld_tableS2.xlxs" |>
    read_langefeld(1)

# GWAS summary statistics
## Langefeld
langefeld_chr2 <- 
    "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Langefeld/ea.imputed.chr2.out" |>
    read_delim(comment = "#", delim = " ") |>
    filter(all_maf >= 0.002) |>
    select(rsid, 
	   pos = position, 
	   alleleA, 
	   alleleB,
	   beta = `frequentist_add_beta_1:add/sle=1`, 
	   se = frequentist_add_se_1, 
	   p = frequentist_add_wald_pvalue_1) |>
    drop_na(beta, se) 

stat4_risk_var <- 
    tier1_ea |>
    filter(gene_region == "STAT4") |>
    slice_min(p) |>
    pull(pos)

langefeld_stat4 <- langefeld_chr2 |>
    filter(between(pos, stat4_risk_var - 5e5, stat4_risk_var + 5e5)) |>
    add_column(chrom = "2", .before = 1) |>
    drop_na(beta, se, p) |>
    mutate(rsid = str_extract(rsid, "rs\\d+")) |>
    select(rsid, pos, alleleA, alleleB, beta, se, p)

# Update and get missing RsIDs by matching position and alleles to dbSNP
dbsnp_meta <- 
    "/lab-share/IM-Gutierrez-e2/Public/References/dbSNP/GCF_000001405.25_GRCh37.p13_assembly_report.txt" |>
    data.table::fread(skip = "# Sequence-Name") 

refseq_chr_id <- dbsnp_meta |>
    filter(`UCSC-style-name` == "chr2") |>
    pull(`RefSeq-Accn`)

dbsnp_stat4_window <- 
    paste0(refseq_chr_id, ":", paste(stat4_risk_var + c(-1e6, 1e6), collapse = "-"))

dbsnp_vcf <- 
    "/lab-share/IM-Gutierrez-e2/Public/References/dbSNP/GCF_000001405.25.gz"

glue("tabix {dbsnp_vcf} {dbsnp_stat4_window}") |>
    paste("| awk '{ print $2,$3,$4,$5 }' > data/dbsnp_stat4.vcf") |>
    system()

dbsnp_stat4 <- "./data/dbsnp_stat4.vcf" |>  
    read_delim(delim = " ", col_names = c("pos", "rsid", "ref", "alt")) |>
    separate_rows(alt, sep = ",")

langefeld_stat4_dbsnp <- 
    langefeld_stat4 |>
    left_join(dbsnp_stat4, join_by(pos, alleleA == ref, alleleB == alt)) |>
    mutate(rsid = ifelse(!is.na(rsid.y), rsid.y, rsid.x)) |>
    select(rsid, pos, alleleA, alleleB, beta, se, p) |>
    filter(!is.na(rsid))

# Save GWAS data
write_tsv(langefeld_stat4_dbsnp, "./data/langefeld_stat4.tsv")

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

