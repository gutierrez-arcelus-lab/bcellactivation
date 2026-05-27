library(tidyverse)
library(susieR)
library(glue)
library(furrr)

run_susie <- function(gene_region, window_idx, summ_stats) {
    
    ld_vcf <- 
	glue("./data/ld/window{window_idx}.vcf.gz") |>
	data.table::fread(skip = "#CHROM") |>
	as_tibble() |>
	unite("ID", c(ID, REF, ALT), sep = "-")

    ld_plink <- 
	glue("./data/ld/window{window_idx}_r.ld") |>
	data.table::fread() |>
	data.matrix()

    rownames(ld_plink) <- colnames(ld_plink) <- ld_vcf$ID
  
    bad_snps <- 
	apply(ld_plink, 2, function(x) all(is.na(x))) |>
	which()

    if ( length(bad_snps) )  ld_plink <- ld_plink[-bad_snps, -bad_snps]

    stats <- 
	summ_stats |>
	filter(locus == gene_region) |>
	select(chrom, pos, rsid, ref, alt, p, beta, se) |>
	mutate(snp_id = glue("{rsid}-{ref}-{alt}")) |>
	filter(snp_id %in% colnames(ld_plink))

    ld_plink_final <- 
	ld_plink[stats$snp_id, stats$snp_id]

    stats_final <- 
	stats |>
	mutate(snp_id = factor(snp_id, levels = colnames(ld_plink))) |>
	arrange(snp_id) |>
	select(chrom, pos, snp_id, p, beta, se) |>
	mutate(z = beta/se)
	
    condz <- 
	kriging_rss(stats_final$z, ld_plink_final, n = sample_size, r_tol = 1e-04)
    
    condz_df <- 
	condz[[2]] |>
	as_tibble() |>
	bind_cols(select(stats_final, snp_id)) |>
	mutate(flag = ifelse(logLR > 2 & abs(z) > 2, TRUE, FALSE))

    write_tsv(condz_df, glue("./data/susie/diagnostics/window{window_idx}.tsv"))

    # if it does not converge, remove problematic SNPs and repeat until convergence
    iter <- 0L
    fit <- list()
    fit$converged <- FALSE
    summ_stats_susie <- stats_final
    ldmat_susie <- ld_plink_final

    while ( !fit$converged & iter <= 30 ) {

	if ( iter > 0 ) { 
	
	    condz <- 
		kriging_rss(summ_stats_susie$z, 
			    ldmat_susie, 
			    n = sample_size,
			    r_tol = 1e-04)

	    snp_rm <- 
		as_tibble(condz$conditional_dist) |>
		rowid_to_column() |>
		arrange(desc(abs(z_std_diff))) |>
		slice(1) |>
		pull(rowid)

	    ldmat_susie <- ldmat_susie[-snp_rm, -snp_rm]
	    summ_stats_susie <- slice(summ_stats_susie, -snp_rm)
	}

	fit <- 
	    susie_rss(summ_stats_susie$z, 
		      ldmat_susie, 
		      n = sample_size, 
		      L = 5, 
		      coverage = 0.9, 
		      r_tol = 1e-05)
    
	write_lines(fit$converged, glue("./data/susie/log/window{window_idx}.txt"), append = TRUE)

	iter <- iter + 1L
    }
    
    write_lines(fit$converged, glue("./data/susie/log/window{window_idx}.txt"))

    if ( !is.null(fit$sets$coverage) ) {

	coverage_df <- 
	    tibble(coverage = scales::percent(round(fit$sets$coverage, 3)),
		   cs = paste0("L", fit$sets$cs_index))

	cs_df <- 
	    enframe(fit$sets$cs, "cs", "rowid") |>
	    unnest(cols = rowid) |>
	    left_join(coverage_df, by = "cs")

	pip_df <- 
	    enframe(fit$pip, "ID", "pip") |>
	    rowid_to_column() |>
	    separate(ID, c("rsid", "ref", "alt"), sep = "-", convert = TRUE) |>
	    left_join(cs_df, join_by(rowid))

    } else {

	pip_df <-
	    enframe(fit$pip, "ID", "pip") |>
	    rowid_to_column() |>
	    separate(ID, c("rsid", "ref", "alt"), sep = "-", convert = TRUE) |>
	    mutate(cs = NA)
    }

    write_tsv(pip_df, glue("./data/susie/pip/window{window_idx}.tsv"))

    fit$lbf_variable |>
	as_tibble() |>
	tibble::rowid_to_column("cs") |>
	mutate(cs = paste0("L", cs)) |>
	pivot_longer(-cs, names_to = "snp_id") |>
	pivot_wider(names_from = cs, values_from = value) |>
	write_tsv(glue("./data/susie/lbf/window{window_idx}.tsv"))
}


# Bentham et al sample size
sample_size <- 6748 + 11516

windows_df <- 
    "./data/windows.tsv" |>
    read_tsv(col_names = c("locus", "window"))

summary_stats <- 
    read_tsv("./data/summary_stats.tsv") |>
    select(locus = gene_region, chrom = chr, rsid, pos, ref = alleleA, alt = alleleB, beta, se, p)

plan(multisession, workers = availableCores())

future_walk2(windows_df$locus, 1:nrow(windows_df), 
	    ~run_susie(gene_region = .x, window_idx = .y, summ_stats = summary_stats),
	    seed = NULL)
