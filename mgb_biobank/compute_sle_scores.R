library(tidyverse)

list.files("./sle_variants")

info_df <- "./sle_variants/sle_variants.tsv" %>%
    read_tsv() %>%
    mutate(or = ifelse(or < 1L, 1L/or, or))

ld_vars <- read_tsv("./sle_variants/sle_ld.tsv") %>%
    filter(r2 >= 0.6) %>%
    distinct(bentham) %>%
    pull(bentham)

var_df <- "./sle_variants/sle_langefeld_bentham_hg38.bed" %>%
    read_tsv(col_names = FALSE) %>%
    separate(X4, c("snp_id", "study"), sep = "-") %>%
    select(chr = 1, pos = 2, snp_id, study) %>%
    left_join(info_df, by = c("chr", "study", "snp_id")) %>%
    filter(!(study == "bentham" & snp_id %in% ld_vars))

vcf <- read_tsv("./sle_variants/sle.MGB.vcf", comment = "##")

vcf_long <- vcf %>%
    select(-(REF:FORMAT)) %>%
    select(chr = 1, pos = POS, snp_id = ID, everything()) %>%
    inner_join(select(var_df, chr, pos)) %>%
    pivot_longer(-(chr:snp_id), names_to = "sample_id", values_to = "genotype")

vcf_long %>% count(genotype)

vcf_recode <- vcf_long %>%
    mutate(dose = case_when(genotype == "0|0" ~ "0",
			    genotype == "1|0" ~ "1",
			    genotype == "0|1" ~ "1",
			    genotype == "1|1" ~ "2",
			    genotype == "0|2" ~ "1",
			    TRUE ~ NA_character_)) %>%
    group_by(snp_id) %>%
    filter(!any(is.na(dose))) %>%
    ungroup()

dosage_df <- vcf_recode %>%
    left_join(select(var_df, chr, pos, or)) %>%
    mutate(beta = log(or))

out <- dosage_df %>%
    group_by(sample_id) %>%
    summarise(het_score = sum(as.integer(dose == "1")),
	      het_score_wt = sum(as.integer(dose == "1") * beta))

write_tsv(out, "./sle_variants/scores.tsv")
