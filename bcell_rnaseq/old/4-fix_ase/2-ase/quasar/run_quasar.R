library(tidyverse)
library(QuASAR)
library(qvalue)

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
    ase_out <- ase_results |>
	map_df("dat", .id = "sample_id") |>
	as_tibble() |>
	select(sample_id, 
	       snp_id = annotations.rsID, 
	       beta = betas, 
	       beta_se = betas.se, 
	       pval = matches("^pval"))
}

get_gt_probs <- function(input_files, min_cov = 10) {

    ase_dat <- UnionExtractFields(input_files, combine = TRUE)

    ase_dat_gt <- PrepForGenotyping(ase_dat, min.coverage = min_cov)

    sample_names <- colnames(ase_dat_gt$ref)

    ase_joint <- 
	fitAseNullMulti(ase_dat_gt$ref, 
			ase_dat_gt$alt, 
			log.gmat = log(ase_dat_gt$gmat))
    
    #allele counts
    ref_counts <-
	bind_cols(ase_dat_gt$ref, ase_dat_gt$annotations) |>
	as_tibble() |>
	pivot_longer(-(chr:af), names_to = "sample_id") |>
	group_by(rsID) |>
	summarise(value = sum(value)) |>
	ungroup() |>
	select(rsID, ref_count = value)
    
    alt_counts <-
	bind_cols(ase_dat_gt$alt, ase_dat_gt$annotations) |>
	as_tibble() |>
	pivot_longer(-(chr:af), names_to = "sample_id") |>
	group_by(rsID) |>
	summarise(value = sum(value)) |>
	ungroup() |>
	select(rsID, alt_count = value)
   
    count_df <- 
	left_join(ref_counts, alt_counts, join_by(rsID))

    #probabilities
    gt_probs <- 
        bind_cols(ase_dat_gt$annotations, ase_joint$gt) |>
        as_tibble() |>
        select(rsID, g0:g2)

    out <- 
	left_join(count_df, gt_probs, join_by(rsID)) |>
	rename("snp_id" = rsID)
    
    out
}


# Pileup files
metadata <- 
    "/lab-share/IM-Gutierrez-e2/Public/vitor/lupus/bcell_rnaseq/4-fix_ase/1-mapping/metadata.tsv" |>
    read_tsv(col_names = c("donor_id", "sample_id", "stim"), col_types = "ccc--") |>
    select(donor_id, sample_id, stim) |>
    mutate(prefix = paste(sample_id, stim, sep = "_")) |>
    select(donor_id, prefix) |>
    {function(x) split(x, x$donor_id)}() |>
    map("prefix")

files_list <- 
    metadata |>
    map(function(x) {
	    files <- sprintf(file.path(temp_dir, "quasar/%s.quasar.in.gz"), x)
	    setNames(files, x)})

# Run quasar
out <- 
    map_df(files_list, run_quasar, .id = "donor_id") |>
    mutate(qval = qvalue(pval)$qvalues)

# Save results
write_tsv(out, "./quasar_results.tsv")


# Genotype probabilities
probs <- map_df(files_list, get_gt_probs, .id = "donor_id")

write_tsv(probs, "./quasar_genotypes.tsv")
