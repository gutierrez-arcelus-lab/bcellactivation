library(tidyverse)

# Allele frequencies among European women in MGB Biobank
snp_af <- 
    "./data/snps_af.bed" |>
    read_tsv(col_names = c("chr", "pos0", "pos", "variant_id", "ref", "alt", "af"))

# Normal-input RNA-seq meta data create for ASEReadCounter analysis
meta <- 
    "../array_spec.tsv" |>
    read_tsv(col_names = c("donor_id", "sample_id", "stim", "mgbid"), 
	     col_types = c(.default = "c")) |>
    unite("id", c("sample_id", "stim"), sep = "_", remove = FALSE) |>
    select(id, donor_id, sample_id, stim)
    
# Convert ASEReadcounter to quasar format
convert <- function(x) {
    sprintf("../results/%s.asereadcounter.txt", x) |>
    read_tsv() |>
    mutate(pos0 = position - 1L) |>
    left_join(snp_af, 
	      join_by(contig == chr, pos0, position == pos, 
		      variantID == variant_id, refAllele == ref, altAllele == alt)) |>
    select(contig, pos0, position, refAllele, altAllele, variantID, af, 
	   refCount, altCount, otherBases)
}

walk(meta$id, 
     ~convert(.x) |> 
     write_tsv(sprintf("./data/%s.quasar.in.gz", .x), col_names = FALSE))
