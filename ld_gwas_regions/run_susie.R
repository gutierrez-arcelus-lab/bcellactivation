library(tidyverse)
library(susieR)
library(glue)
library(furrr)

run_susie <- function(gene, region_coords, summ_stats) {
    
    ld_vcf <- 
	glue("./data/chr{region_coords}.vcf.gz") |>
	data.table::fread(skip = "#CHROM") |>
	as_tibble() |>
	unite("ID", c(ID, REF, ALT), sep = "-")

    ld_plink <- 
	glue("./data/chr{region_coords}_r.ld") |>
	data.table::fread() |>
	data.matrix()

    rownames(ld_plink) <- colnames(ld_plink) <- ld_vcf$ID

    stats <- 
	summ_stats |>
	filter(gene_region == gene) |>
	select(chr, pos, rsid, alleleA, alleleB, p, beta, se) |>
	mutate(snp_id = glue("{rsid}-{alleleA}-{alleleB}")) |>
	filter(snp_id %in% colnames(ld_plink))

    ld_plink_final <- ld_plink[stats$snp_id, stats$snp_id]

    stats_final <- 
		stats |>
		mutate(snp_id = factor(snp_id, levels = colnames(ld_plink))) |>
		arrange(snp_id) |>
		select(chr, pos, snp_id, p, beta, se) |>
		mutate(z = beta/se)
	
    condz <- kriging_rss(stats_final$z, ld_plink_final, n = sample_size, r_tol = 1e-04)
    
    condz_df <- 
		condz[[2]] |>
		as_tibble() |>
		bind_cols(select(stats_final, snp_id)) |>
		mutate(flag = ifelse(logLR > 2 & abs(z) > 2, TRUE, FALSE))

    write_tsv(condz_df, glue("./data/susie_{gene}_krss.tsv"))

    # if it does not converge, remove problematic SNPs and repeat until convergence
    iter <- 0L
    fit <- list()
    fit$converged <- FALSE
    summ_stats_susie <- stats_final
    ldmat_susie <- ld_plink_final

    while ( !fit$converged & iter <= 30 ) {

	if ( iter > 0 ) { 
	
	    condz <- kriging_rss(summ_stats_susie$z, 
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
			 L = 3, 
			 coverage = 0.9, 
			 r_tol = 1e-05)

		iter <- iter + 1L
    }
    
    write_lines(fit$converged, glue("./data/susie_{gene}_converged.txt"))

    if (! is.null(fit$sets$coverage) ) {

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
	    left_join(cs_df, join_by(rowid)) |>
	    mutate(cs = ifelse(!is.na(cs), paste0(cs, " (", coverage, ")"), NA))
    } else {

	pip_df <-
	    enframe(fit$pip, "ID", "pip") |>
	    rowid_to_column() |>
	    separate(ID, c("rsid", "ref", "alt"), sep = "-", convert = TRUE) |>
	    mutate(cs = NA)
    }

    write_tsv(pip_df, glue("./data/susie_{gene}_pip.tsv"))
}


# Langefeld et al sample size
sample_size <- 4036 + 6959

regions <- "./data/langefeld_regions.tsv" |>
    read_tsv(col_names = c("region", "gene"))

summary_stats <- "./data/langefeld_summ_stats.tsv" |>
    read_tsv()


plan(multisession, workers = availableCores())

future_walk2(regions$gene, regions$region, 
	    ~run_susie(gene = .x, region_coords = .y, summ_stats = summary_stats),
	    seed = NULL)

