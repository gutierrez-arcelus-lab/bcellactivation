library(tidyverse)

bed_hg19 <- 
    read_tsv("./sle_variants/sle.bed", col_names = FALSE) %>%
    select(chr = X1, pos_hg19 = X3, snp_id = X4)

bed_hg38 <- 
    read_tsv("./sle_variants/sle_hg38.bed", col_names = FALSE) %>%
    select(chr = X1, pos_hg38 = X3, snp_id = X4)

bed <- left_join(bed_hg19, bed_hg38, by = c("chr", "snp_id")) %>%
    select(-snp_id)

# If OR < 0, use the reciprocal for risk
sle_vars <- read_tsv("./sle_variants/sle.tsv") %>%
    select(pos, or) %>%
    mutate(or = ifelse(or < 1L, 1L/or, or))

# import VCF

vcf <- read_tsv("./sle_variants/sle.MGB.vcf", comment = "##") %>%
    filter(nchar(REF) == 1L, nchar(ALT) == 1L) %>%
    select(-(QUAL:FORMAT)) %>%
    select(chr = 1, pos = POS, snp_id = ID, ref = REF, alt = ALT, everything()) %>%
    pivot_longer(-(chr:alt), names_to = "sample_id", values_to = "genotype") %>%
    mutate(dose = case_when(genotype == "0|0" ~ "0",
			    genotype == "1|0" ~ "1",
			    genotype == "0|1" ~ "1",
			    genotype == "1|1" ~ "2",
			    TRUE ~ NA_character_))

dosage_df <- vcf %>%
    left_join(bed, by = c("chr", "pos" = "pos_hg38")) %>%
    left_join(sle_vars, by = c("pos_hg19" = "pos")) %>%
    mutate(beta = log(or))

out <- dosage_df %>%
    group_by(sample_id) %>%
    summarise(het_score = sum(as.integer(dose == "1")),
	      het_score_wt = sum(as.integer(dose == "1") * beta))

write_tsv(out, "./sle_variants/scores.tsv")
