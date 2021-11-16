library(tidyverse)


bed_hg19 <- 
    read_tsv("./sle_variants/sle.bed", col_names = FALSE) %>%
    select(chr = X1, pos_hg19 = X3, snp_id = X4)

bed_hg38 <- 
    read_tsv("./sle_variants/sle_hg38.bed", col_names = FALSE) %>%
    select(chr = X1, pos_hg38 = X3, snp_id = X4)

bed <- left_join(bed_hg19, bed_hg38, by = c("chr", "snp_id")) %>%
    select(-snp_id)

sle_vars <- read_tsv("./sle_variants/sle.tsv") %>%
    select(pos, or)

vcf <- read_tsv("./sle_variants/sle.MGB.vcf", comment = "##") %>%
    select(-(QUAL:FORMAT)) %>%
    select(chr = 1, pos = POS, snp_id = ID, ref = REF, alt = ALT, everything()) %>%
    pivot_longer(-(chr:alt), names_to = "sample_id", values_to = "genotype") %>%
    separate(genotype, c("h1", "h2"), sep = "\\|") %>%
    pivot_longer(h1:h2, names_to = "haplot", values_to = "allele") %>%
    mutate(allele = as.integer(allele))

dosage_df <- vcf %>%
    left_join(bed, by = c("chr", "pos" = "pos_hg38")) %>%
    left_join(sle_vars, by = c("pos_hg19" = "pos")) %>%
    mutate(beta = abs(log(or))) %>%
    group_by(chr, pos, snp_id, ref, alt, sample_id, beta) %>%
    summarise(dose = sum(allele)) %>%
    ungroup()

out <- dosage_df %>%
    group_by(sample_id) %>%
    summarise(het_score = sum(as.integer(dose == 1L)),
	      het_score_wt = sum(as.integer(dose == 1) * beta)) %>%
    ungroup()

write_tsv(out, "./sle_variants/scores.tsv")
