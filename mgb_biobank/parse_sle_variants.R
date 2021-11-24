library(tidyverse)
library(readxl)

read_tab <- function(tab_path, sheet_number) {
    read_excel(tab_path, sheet_number, skip = 2) %>%
    select(snp_id = 1, 
	   chr = 2, 
	   pos = 3, 
	   gene_region = 4, 
	   region_rank = 5, 
	   ref_allele = 6, 
	   p = `P-value`,
	   or = `OR (95% CI)`,
	   p_stepwise = `Regional Stepwise  P-value`) %>% 
    filter(!is.na(pos)) %>%
    mutate(p = gsub("\\s+", "", p), 
	   p = sub("(10-)(\\d+)", "1e-\\2", p),
	   or = parse_number(or)) %>%
    separate(p, c("n", "p"), sep = "x") %>%
    mutate(n = as.numeric(n),
	   p = as.numeric(p),
	   p_value = n * p) %>%
    select(snp_id, chr, pos, gene_region, region_rank, ref_allele, p_value, p_stepwise, or)
}

tier1_ea <- read_tab("./sle_variants/langefeld_tableS2.xlxs", 1)
tier2_ea <- read_tab("./sle_variants/langefeld_tableS2.xlxs", 2)
tier3_ea <- read_tab("./sle_variants/langefeld_tableS2.xlxs", 3)

tab <- bind_rows("1" = tier1_ea, "2" = tier2_ea, "3" = tier3_ea, .id = "tier") %>%
    arrange(tier, region_rank, p_value)

tab_best <- tab %>% 
    group_by(tier, region_rank) %>% 
    slice(which.min(p_value)) %>%
    ungroup()

#tab_fdr01 <- tab %>%
#    filter(p_value <= 1/46744 * 0.01)

bed <- tab_best %>%
    mutate(chrom = sub("^(\\d+).*$", "chr\\1", chr),
	   chrom = factor(chrom, levels = paste0("chr", 1:22)),
	   snp_id = sub("^(\\S+)\\s+.*$", "\\1", snp_id)) %>%
    mutate(start = pos - 1L) %>% 
    select(chrom, start, end = pos, snp_id) %>%
    arrange(chrom, start) 

#write_tsv(tab_fdr01, "./sle_variants/sle.tsv")
write_tsv(tab_best, "./sle_variants/sle.tsv")
write_tsv(bed, "./sle_variants/sle.bed", col_names = FALSE)
