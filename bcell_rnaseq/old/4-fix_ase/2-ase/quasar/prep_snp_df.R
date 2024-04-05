library(tidyverse)

snp_info <- 
    "/temp_work/ch229163/vcf/quasar/chr%s.merged.info.vcf.gz" |>
    sprintf(c(1:22, "X")) |>
    map_dfr(~read_tsv(., comment = "##") |>
	    mutate(af = str_extract(INFO, "(?<=AF=)[^;]+"),
		   hwe = str_extract(INFO, "(?<=HWE=)[^;]+")) |>
	    mutate_at(vars(af, hwe), as.numeric) |>
	    select(chr = `#CHROM`, pos = POS, snp_id = ID, ref = REF, alt = ALT, af, hwe))

bed <- 
    snp_info |>
    filter(hwe > 0.05/nrow(snp_info)) |>
    mutate(start = pos - 1L) |>
    select(chr, start, pos, snp_id, ref, alt, af)

write_tsv(bed, "./data/snps_af.bed", col_names = FALSE)
