library(tidyverse)

pileup <- read_tsv("./demuxlet/pileup.txt", col_names = FALSE)

pile_counts <- pileup %>%
    mutate(matchbases = str_extract_all(X5, "[AGCTagct]")) %>%
    select(chr = X1, pos = X2, matchbases) %>%
    unnest(cols = c(matchbases)) %>%
    mutate(matchbases = toupper(matchbases)) %>%
    count(chr, pos, matchbases)

vcf <- read_tsv("./demuxlet/samples_allsnps.vcf.gz", comment = "##") %>%
    select(chr = 1, pos = POS, ID, REF, ALT, SAMPLE1, SAMPLE2)

match_df <- vcf %>%
    inner_join(pile_counts) %>%
    separate(SAMPLE1, c("SAMPLE1.1", "SAMPLE1.2"), sep = "/") %>%
    separate(SAMPLE2, c("SAMPLE2.1", "SAMPLE2.2"), sep = "/") %>%
    mutate_at(vars(SAMPLE1.1:SAMPLE2.2), ~ifelse(. == 0, REF, ALT)) %>%
    select(chr, pos, ID, SAMPLE1.1:SAMPLE2.2, matchbases) %>%
    mutate_at(vars(SAMPLE1.1:SAMPLE2.2), ~ifelse(. == matchbases, 1L, 0L)) %>%
    select(-matchbases) %>%
    pivot_longer(SAMPLE1.1:SAMPLE2.2, names_to = "hap", values_to = "match") %>%
    separate(hap, c("sample", "hap"), sep = "\\.")
    
mismatch_snps <- match_df %>%
    group_by(chr, pos, ID) %>%
    filter(all(match == 0)) %>%
    ungroup() %>%
    distinct(chr, pos, ID)

mismatch_snps %>%
    inner_join(pile_counts) %>%
    summarise(total = sum(n))


