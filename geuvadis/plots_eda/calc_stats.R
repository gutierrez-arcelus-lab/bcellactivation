library(tidyverse)

return_vcf_types <- function(vcf_in)     
    vcf_in %>%
	rowid_to_column() %>%
	separate_rows(alt, sep = ",") %>%
	mutate(var_len = pmax(nchar(ref), nchar(alt))) %>%
	group_by(across(c(-alt))) %>%
	summarise(len = max(var_len),
		  alt = paste(alt, collapse = ",")) %>%
	ungroup() %>%
	mutate(type = ifelse(len > 1L, "INDEL", "SNV")) %>%
	select(-rowid, -len, -var_len)

vcf <- read_tsv("./ERR188022.norm.vcf.gz", comment = "##") %>%
    separate(ERR188022, c("gt", "ad", "dp", "gq", "pl"), sep = ":") %>%
    select(chr = 1, pos = 2, ref = 4, alt = 5, gt, dp) %>%
    filter(ref != "*" & alt != "*")

vcf_types <- return_vcf_types(vcf)

summ_types <- vcf_types %>%
    count(type) %>%
    mutate(p = n/sum(n) * 100) %>%
    select(class = 1, n, p)

summ_snv_class <- vcf_types %>%
    filter(type == "SNV") %>%
    mutate(biallelic = ifelse(!grepl(",", alt), "biallelic", "multiallelic")) %>%
    count(biallelic) %>%
    mutate(p = round(n/sum(n) * 100, 2)) %>%
    select(class = 1, n, p)

summ_gt <- vcf_types %>%
    filter(type == "SNV") %>%
    count(gt, sort = TRUE) %>%
    mutate(p = round(n/sum(n) * 100, 2)) %>%
    select(class = 1, n, p)

bind_rows("variant type" = summ_types,
	  "SNV class" = summ_snv_class,
	  "Biallelic SNV genotype" = summ_gt,
	  .id = "stat") %>%
    write_tsv("summary_stats.tsv")

vcf %>%
    pull(dp) %>%
    write_lines("./depth.txt")


# Comparison to 1000G data
vcf_1000g_lowcov <- "../results/gatk/ERR188022_1000KG.vcf.gz" %>%
    read_tsv(comment = "##") %>%
    select(chr = 1, pos = 2, ref = 4, alt = 5, gt = NA12812) %>%
    mutate(chr = paste0("chr", chr)) %>%
    return_vcf_types()

isec_1000g <- left_join(vcf_types, vcf_1000g_lowcov, 
	  by = c("chr", "pos", "type"),
	  suffix = c("_rna", "_low1000g"))

summ_1000g_types <- count(isec_1000g, type, present = !is.na(ref_low1000g)) %>%
    mutate(present = recode(as.character(present), "TRUE" = "present", "FALSE" = "NA")) %>%
    group_by(type) %>%
    mutate(p = n/sum(n) * 100) %>%
    ungroup() %>%
    unite("category", c("type", "present"), sep = "_")

isec_1000g_snp <- isec_1000g %>%
    filter(type == "SNV", !is.na(ref_low1000g), !grepl(",", alt_rna)) %>%
    mutate(gt_rna = ifelse(gt_rna == "1/0", "0/1", gt_rna),
	   gt_low1000g = sub("\\|", "/", gt_low1000g),
	   gt_low1000g = ifelse(gt_low1000g == "1/0", "0/1", gt_low1000g))

summ_1000g_alleles <- isec_1000g_snp %>%
    mutate(category = case_when(ref_rna == ref_low1000g & alt_rna == alt_low1000g ~ "ref & alt",
				ref_rna == ref_low1000g & alt_rna != alt_low1000g ~ "diff alt",
				TRUE ~ "other")) %>%
    count(category) %>%
    mutate(p = n/sum(n) * 100)

summ_1000g_gts <- isec_1000g_snp %>%
    count(gt_rna, gt_low1000g, sort = TRUE) %>%
    unite("category", c("gt_rna", "gt_low1000g"), sep = " - ") %>%
    mutate(p = n/sum(n) * 100)

summ_1000g <- bind_rows("variant type" = summ_1000g_types,
			"bi-SNV REF/ALT match" = summ_1000g_alleles,
			"bi-SNV Genotypes" = summ_1000g_gts, 
			.id = "stat")

# Comparison to 1000G high coverage data
vcf_1000g_high <- "../results/gatk/ERR188022.nygc1000g.vcf.gz" %>%
    read_tsv(comment = "##") %>%
    select(chr = 1, pos = 2, ref = 4, alt = 5, gt = NA12812) %>%
    return_vcf_types()

isec_1000g_high <- left_join(vcf_types, vcf_1000g_high, 
	  by = c("chr", "pos", "type"),
	  suffix = c("_rna", "_high1000g"))

summ_high1000g_types <- count(isec_1000g_high, type, present = !is.na(ref_high1000g)) %>%
    mutate(present = recode(as.character(present), "TRUE" = "present", "FALSE" = "NA")) %>%
    group_by(type) %>%
    mutate(p = n/sum(n) * 100) %>%
    ungroup() %>%
    unite("category", c("type", "present"), sep = "_")

isec_high1000g_snp <- isec_1000g_high %>%
    filter(type == "SNV", !is.na(ref_high1000g), !grepl(",", alt_rna)) %>%
    mutate(gt_rna = ifelse(gt_rna == "1/0", "0/1", gt_rna),
	   gt_high1000g = sub("\\|", "/", gt_high1000g),
	   gt_high1000g = ifelse(gt_high1000g == "1/0", "0/1", gt_high1000g))

summ_high1000g_alleles <- isec_high1000g_snp %>%
    mutate(category = case_when(ref_rna == ref_high1000g & alt_rna == alt_high1000g ~ "ref & alt",
				ref_rna == ref_high1000g & alt_rna != alt_high1000g ~ "diff alt",
				TRUE ~ "other")) %>%
    count(category) %>%
    mutate(p = n/sum(n) * 100)

summ_high1000g_gts <- isec_high1000g_snp %>%
    count(gt_rna, gt_high1000g, sort = TRUE) %>%
    unite("category", c("gt_rna", "gt_high1000g"), sep = " - ") %>%
    mutate(p = n/sum(n) * 100)

summ_1000g_high <- bind_rows("variant type" = summ_high1000g_types,
			     "bi-SNV REF/ALT match" = summ_high1000g_alleles,
			     "bi-SNV Genotypes" = summ_high1000g_gts, 
			     .id = "stat")

left_join(summ_1000g, summ_1000g_high, 
	  by = c("stat", "category"),
	  suffix = c("_lowcov", "_highcov")) %>%
    write_tsv("summary_1000g_match.tsv")


