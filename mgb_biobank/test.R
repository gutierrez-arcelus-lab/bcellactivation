library(tidyverse)


vcf <- read_tsv("./sle_variants/sle.MGB.vcf", comment = "##") %>%
    select(-(QUAL:FORMAT)) %>%
    select(chr = 1, pos = POS, snp_id = ID, ref = REF, alt = ALT, everything()) %>%
    pivot_longer(-(chr:alt), names_to = "sample_id", values_to = "genotype") %>%
    separate(genotype, c("h1", "h2"), sep = "\\|") %>%
    pivot_longer(h1:h2, names_to = "haplot", values_to = "allele") %>%
    mutate(allele = as.integer(allele))

dosage_df <- vcf %>%
    group_by(chr, pos, snp_id, ref, alt, sample_id) %>%
    summarise(dose = sum(allele)) %>%
    ungroup()

dosage_df %>%
    group_by(sample_id) %>%
    summarise(het_score = sum(as.integer(dose == 1L))) %>%
    ungroup() %>%
    arrange(desc(het_score)) %>% 
    top_n(12, het_score)
