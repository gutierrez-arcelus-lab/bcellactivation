library(tidyverse)
library(readxl)

tab <- "./sle_variants/langefeld_tableS2.xlxs" %>%
    read_excel(3, skip = 2) %>%
    select(snp_id = 1, chr = 2, pos = 3, gene = 4, 
	   region_rank = 5, ref_allele = 6, p = `P-value`,
	   p_stepwise = `Regional Stepwise  P-value`) %>% 
    filter(!is.na(pos)) %>%
    mutate(p = gsub("\\s+", "", p), 
	   p = sub("(10-)(\\d+)", "1e-\\2", p)) %>%
    separate(p, c("n", "p"), sep = "x") %>%
    mutate(n = as.numeric(n),
	   p = as.numeric(p),
	   p_value = n * p) %>%
    select(snp_id, chr, pos, gene, region_rank, ref_allele, p_value, p_stepwise)

