library(Biostrings)
library(tidyverse)

fasta <- readDNAStringSet("./GRCh38.primary_assembly.genome.fa")

names(fasta) <- sub(" (.+)$", "", names(fasta))

writeXStringSet(fasta, "./GRCh38.primary_assembly.genome.fixnames.fa")

gwas <- read_tsv("./gwas_catalog_v1.0.2.tsv", guess_max = 1e5)

bentham <- gwas %>%
    filter(grepl("lupus", `DISEASE/TRAIT`, ignore.case = TRUE)) %>%
    filter(grepl("Bentham", `FIRST AUTHOR`)) %>%
    select(`FIRST AUTHOR`, CHR_ID, CHR_POS, `REPORTED GENE(S)`, SNPS, SNP_ID_CURRENT,
	   `P-VALUE`, `PVALUE_MLOG`, `OR or BETA`) %>%
    mutate(CHR_ID = factor(CHR_ID, levels = c(1:22, "X"))) %>%
    arrange(`FIRST AUTHOR`, CHR_ID, CHR_POS)

bentham %>%
    group_by(`REPORTED GENE(S)`) %>%
    dplyr::slice(which.min(`P-VALUE`)) %>%
    ungroup() %>%
    arrange(CHR_ID, CHR_POS) 

gwas %>% 
    mutate(CHR_POS = as.numeric(CHR_POS)) %>%
    filter(CHR_ID == 1, between(CHR_POS, 207479600, 207479900)) %>%
    select(`FIRST AUTHOR`, CHR_ID, CHR_POS, `REPORTED GENE(S)`, SNPS)

gwas %>%
    filter(grepl("rs2182909", SNPS))

grep("gene", names(gwas), ignore.case = T, value=T)

gwas %>%
    select(`FIRST AUTHOR`, CHR_ID, CHR_POS, SNPS, contains("GENE"), contains("TRAIT")) %>%
    separate_rows(`REPORTED GENE(S)`, sep = ",") %>%
    mutate(`REPORTED GENE(S)` = trimws(`REPORTED GENE(S)`)) %>%
    filter(`REPORTED GENE(S)` == "CR2") %>%
    as.data.frame()





##


gwas |>
    filter(grepl("Ferreira", `FIRST AUTHOR`)) |> 
    filter(`DISEASE/TRAIT` == "Asthma (adult onset)") |>
    print(width = Inf, n = Inf)
