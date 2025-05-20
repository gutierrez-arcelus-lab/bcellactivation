library(tidyverse)
library(furrr)

snp_info <- 
    read_tsv("./data/Khunsriraksakul_summstats_hg38.tsv") |>
    select(chr, pos, rsid, ref, alt)

sentinels <- read_tsv("./data/regions.tsv", col_names = c("region", "rsid"))

temp_dir <- system("echo $TEMP_WORK", intern = T)

main <- function(s_snp) {
    
    vcf <- 
	file.path(temp_dir, "vcf/prs/%s.filt.vcf.gz") |>
	sprintf(s_snp) |>
	read_tsv(comment = "##") 

    plink <- 
	file.path(temp_dir, "vcf/prs/%s.ld") |>
	sprintf(s_snp) |>
	read_delim(delim = " ", col_names = FALSE)

    if ( (ncol(plink) %% nrow(plink) == 1) && all(is.na(plink[, ncol(plink)])) ) {

	plink <- plink[, -ncol(plink)]
    }

    valid_snps <- which(sapply(plink, function(x) !all(is.na(x))))

    plink_filt <- plink[valid_snps, valid_snps]
    ldmat <- data.matrix(plink_filt)
    vcf <- slice(vcf, valid_snps)
    rownames(ldmat) <- vcf$ID
    colnames(ldmat) <- vcf$ID

    ldmat_filt <- 
	ldmat |>
	as_tibble(rownames = "id1") |>
	pivot_longer(-id1, names_to = "id2", values_to = "r2") |>
	mutate(idmin = pmin(id1, id2),
	       idmax = pmax(id1, id2)) |>
	distinct(idmin, idmax, r2) |>
	filter(idmin != idmax, r2 >= .8)

    if ( nrow(ldmat_filt) == 0 ) {
	return(tibble(sentinel = s_snp, snp2 = NA, r2 = NA)) 
    }

    # Do inner_join at the end to filter for SNPs included in summary stats
    # because filter with bcftools also include indels that overlap a position
    tmp <- 
	ldmat_filt |>
	select(idmin, idmax) |>
	pivot_longer(idmin:idmax) |>
	distinct(snp = value) |>
	separate(snp, 
		 c("chr", "pos", "REF", "ALT"), 
		 sep = ":", remove = FALSE, 
		 convert = TRUE) |>
	mutate(chr = paste0("chr", chr)) |>
	inner_join(snp_info, join_by(chr, pos))

    tmp_match_alleles <- 
	tmp |>
	separate_rows(alt, sep = ",") |>
	group_by(snp) |>
	filter(REF == ref, ALT == alt) |>
	ungroup()

    out <- ldmat_filt |>
	left_join(select(tmp_match_alleles, snp, rsid), join_by(idmin == snp)) |>
	select(snp1 = rsid, idmax, r2) |>
	left_join(select(tmp_match_alleles, snp, rsid), join_by(idmax == snp)) |>
	select(snp1, snp2 = rsid, r2) |>
	filter(snp1 == s_snp | snp2 == s_snp) |>
	mutate(sentinel = ifelse(snp1 == s_snp, snp1, snp2),
	       snp2 = ifelse(snp2 == s_snp, snp1, snp2)) |>
	select(sentinel, snp2, r2) |>
	arrange(desc(r2))

    if ( nrow(out) == 0 ) {
	return(tibble(sentinel = s_snp, snp2 = NA, r2 = NA)) 
    } else {
	return(out)
    }
}

# run analysis
plan(multisession, workers = availableCores())

res <- future_map(sentinels$rsid, safely(~main(s_snp = .x)))

errors <- map(res, "error")

errors |>
	map(~!is.null(.)) |>
	unlist() |>
	which()

res_df <- res |>
    map("result") |>
    bind_rows()

write_tsv(res_df, "./data/Khunsriraksakul_sentinel_r2.tsv")
