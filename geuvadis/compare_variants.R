library(tidyverse)

get_vartype <- function(REF, ALT) {

    ref_indel <- str_split(REF, ",") %>% unlist() %>% nchar %>% `>`(1)
    alt_indel <- str_split(ALT, ",") %>% unlist() %>% nchar %>% `>`(1) %>% any()

    case_when((ref_indel == TRUE || alt_indel == TRUE) ~ "INDEL",
	      TRUE ~ "SNP")
}


rna_vcf <- "./results/gatk/ERR188022.filtered.vcf.gz" %>%
    read_tsv(comment = "##") %>%
    select(chr = 1, pos = 2, ref = 4, alt = 5, filter = 7, geno = ERR188022) %>%
    mutate(geno = sub("^([^:]+).*$", "\\1", geno),
	   var_type = map2_chr(ref, alt, get_vartype))

kgp_vcf <- "./results/gatk/ERR188022_1000KG.vcf.gz" %>%
    read_tsv(comment = "##") %>%
    select(chr = 1, pos = 2, ref = 4, alt = 5, geno = NA12812) %>%
    mutate(chr = paste0("chr", chr),
	   geno = sub("\\|", "/", geno),
	   geno = ifelse(geno == "1/0", "0/1", geno),
	   var_type = map2_chr(ref, alt, get_vartype))

dbsnp_vcf <- "./results/gatk/ERR188022_dbsnp.vcf.gz" %>%
    read_tsv(comment = "##") %>%
    select(chr = 1, pos = 2, ref = 4, alt = 5) %>%
    mutate(var_type = map2_chr(ref, alt, get_vartype))



full_join(rna_vcf, dbsnp_vcf, by = c("chr", "pos", "var_type"),
	  suffix = c("_rna", "_dbsnp")) %>%
left_join(kgp_vcf, by = c("chr", "pos", "var_type"), 
	  suffix = c("_rna", "_1000g")) %>%
mutate(mismatch = case_when(is.na(ref_dbsnp) & is.na(ref) ~ "NA in dbSNP & 1000G",
			    !is.na(ref_dbsnp) & is.na(ref) ~ "NA in 1000G",
			    is.na(ref_dbsnp) & !is.na(ref) ~ "NA in dbSNP",
			    is.na(ref_rna) ~ "NA in RNA-seq call",
			    (ref_rna == ref & alt_rna == alt) ~ "correct",
			    ref_rna == ref_dbsnp & ref_rna == ref & alt_rna != alt ~ "wrong ALT",
			    ref_rna == alt & alt_rna == ref ~ "swap REF-ALT",
			    TRUE ~ "other")) %>%
filter(mismatch == "other") %>%
print(n = 20)

rna_vcf %>%
    filter(chr == "chr1", pos == 84929933)

kgp_vcf %>%
    filter(chr == "chr1", pos == 84929933)

dbsnp_vcf %>%
    filter(chr == "chr1", pos == 84929933)


# de-normalize VCF from RNA-seq
# get variant type from INFO field for dbSNP
# convert 0,1,2 to nucleotides for comparison








