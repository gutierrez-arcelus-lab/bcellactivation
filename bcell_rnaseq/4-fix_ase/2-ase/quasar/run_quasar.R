library(tidyverse)
library(QuASAR)

system("export BEDTOOLS_X=2.31.0")
temp_dir <- system("echo $TEMP_WORK", intern = TRUE) 

# Function
run_quasar <- function(input_files, min_cov = 10) {

    ase_dat <- UnionExtractFields(input_files, combine = TRUE)

    ase_dat_gt <- PrepForGenotyping(ase_dat, min.coverage = min_cov)

    sample_names <- colnames(ase_dat_gt$ref)

    ase_joint <- 
	fitAseNullMulti(ase_dat_gt$ref, 
			ase_dat_gt$alt, 
			log.gmat = log(ase_dat_gt$gmat))

    #gt_probs <- 
    #    bind_cols(ase_dat_gt$annotations, ase_joint$gt) |>
    #    as_tibble() |>
    #    select(snp_id = rsID, g0:g2)

    # Estimate ASE
    ase_results <- 
	aseInference(gts = ase_joint$gt, 
		     eps.vect = ase_joint$eps, 
		     priors = ase_dat_gt$gmat, 
		     ref.mat = ase_dat_gt$ref, 
		     alt.mat = ase_dat_gt$alt, 
		     min.cov = min_cov, 
		     sample.names = sample_names, 
		     annos = ase_dat_gt$annotations)

    # Output
    ase_results |>
	map_df("dat", .id = "sample_id") |>
	select(sample_id, 
	       snp_id = annotations.rsID, 
	       beta = betas, 
	       beta_se = betas.se, 
	       pval = matches("^pval")) |>
	#left_join(gt_probs, join_by(snp_id)) |>
	#select(sample_id, snp_id, g0:g2, beta:pval) |>
	as_tibble()
}

# Pileup files
metadata <- 
    "/lab-share/IM-Gutierrez-e2/Public/vitor/lupus/bcell_rnaseq/4-fix_ase/1-mapping/metadata.tsv" |>
    read_tsv() |>
    select(donor_id, sample_id, stim) |>
    mutate(prefix = paste(sample_id, stim, sep = "_")) |>
    select(donor_id, prefix) |>
    {function(x) split(x, x$donor_id)}() |>
    map("prefix")

files_list <- 
    metadata |>
    map(function(x) {
	    files <- sprintf(file.path(temp, "quasar/%s.quasar.in.gz"), x)
	    setNames(files, x)})

# Run quasar
out <- map_df(files_list, run_quasar, .id = "donor_id")

# Save results
write_tsv(out, "./quasar_results.tsv")
