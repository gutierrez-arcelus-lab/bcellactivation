library(tidyverse)

get_vartype <- function(REF, ALT) {

    ref_indel <- str_split(REF, ",") %>% unlist() %>% nchar %>% `>`(1)
    alt_indel <- str_split(ALT, ",") %>% unlist() %>% nchar %>% `>`(1) %>% any()

    case_when((ref_indel == TRUE || alt_indel == TRUE) ~ "INDEL",
	      TRUE ~ "SNV")
}

prefix <- commandArgs(TRUE)[1]
sampleid <- basename(prefix)

rna_vcf <- paste0(prefix, ".filtered.denorm.vcf.gz") %>%
    read_tsv(comment = "##") %>%
    separate({{ sampleid }}, c("geno", "AD", "depth", "GQ", "PL"), sep = ":", convert = TRUE) %>%
    select(chr = 1, pos = 2, ref = 4, alt = 5, filter = 7, geno, depth) %>%
    mutate(var_type = map2_chr(ref, alt, get_vartype))

#kgp_vcf <- "./results/gatk/ERR188022_1000KG.vcf.gz" %>%
#    read_tsv(comment = "##") %>%
#    select(chr = 1, pos = 2, ref = 4, alt = 5, geno = NA12812) %>%
#    mutate(chr = paste0("chr", chr),
#	   geno = sub("\\|", "/", geno),
#	   geno = ifelse(geno == "1/0", "0/1", geno),
#	   var_type = map2_chr(ref, alt, get_vartype))
#
dbsnp_vcf <- paste0(prefix, "_dbsnp.vcf.gz") %>% 
    read_tsv(comment = "##") %>%
    mutate(var_type = str_extract(INFO, "(?<=VC=)[^;]+")) %>%
    select(chr = 1, pos = 2, ref = 4, alt = 5, var_type)

chr_order <- rna_vcf %>%
    distinct(chr) %>%
    pull(chr) %>%
    str_sort(numeric = TRUE)

filtered_positions <- rna_vcf %>%
    filter(var_type == "SNV", geno == "0/1") %>%
    inner_join(dbsnp_vcf, by = c("chr", "pos", "var_type")) %>%
    select(chr, pos) %>%
    mutate(chr = factor(chr, levels = chr_order)) %>%
    arrange(chr, pos) %>%
    distinct(chr, pos)

write_tsv(filtered_positions, 
	  paste0(prefix, ".mergeDBSNP.pos"),
	  col_names = FALSE)
