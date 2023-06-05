library(tidyverse)

### First MBV analysis
metadata <- read_tsv("../../data/metadata.tsv", col_types = c(.default = "c")) |>
    select(batch, donor_vcf = subject_id, vcfid = mgb_id)

mbv_files <- list.files("../../results/mbv/")

mbv_df <- file.path("../../results/mbv", mbv_files) |>
    setNames(mbv_files) |>
    map_dfr(~read_delim(., delim = " "), .id = "bam") |>
    extract(bam, c("donor_id_bam", "bam_stim"), "(\\d+\\.\\d)_([^.]+)\\.txt") |>
    separate(donor_id_bam, c("donor_bam", "rep_bam"), sep = "\\.", remove = FALSE) |>
    extract(SampleID, c("donor_vcf"), "\\d+_[^-]+-(\\d+)")

to_fix_1 <- mbv_df |>
    group_by(donor_id_bam, bam_stim) |>
    filter(perc_het_consistent == max(perc_het_consistent)) |>
    ungroup() |>
    filter(donor_bam != donor_vcf) |>
    filter(perc_het_consistent > 0.9) |>
    left_join(metadata, join_by(donor_vcf)) |>
    select(donor_bam, bam_stim, donor_vcf, vcfid)


### Second MBV analysis (whole MGB biobank)
meta <- read_tsv("../mbv.spec", col_names = FALSE)

files <- file.path("../results", paste0(meta$X1, "_", meta$X2, ".txt"))
names(files) <- meta$X1

res <- files |> 
    map_dfr(~read_delim(., delim = " "), .id = "bamid") |>
    select(bamid, vcfid = SampleID, perc_het_consistent, perc_hom_consistent)

matches_df <- res |> 
    group_by(bamid) |>
    filter(perc_het_consistent == max(perc_het_consistent)) |>
    ungroup()

to_fix_2 <- matches_df |> 
    separate(bamid, c("donor_id", "stim"), sep = "_") |>
    separate(donor_id, c("donor_id", "rep"), sep = "\\.") |>
    extract(vcfid, "vcf_donor_id", ".+-(\\d+)", remove = FALSE) |>
    filter(donor_id != vcf_donor_id) |>
    select(donor_bam = donor_id, bam_stim = stim, donor_vcf = vcf_donor_id, vcfid)

## Set up fix analysis
to_fix_df <- bind_rows(to_fix_1, to_fix_2) |>
    distinct()

mgb_files <- "/lab-share/IM-Gutierrez-e2/Public/External_datasets/mgb_biobank/04%02d/etc/ids.tsv" |>
    sprintf(c(1:8, 10))

names(mgb_files) <- sprintf("04%02d", c(1:8, 10)) 

batch_ids <- 
    mgb_files |>
    map_dfr(~read_tsv(., col_names = c("donor_id", "vcf_id")),
	    .id = "batch")

batch_ids |> 
    filter(vcf_id %in% unique(to_fix_df$vcfid)) |>
    distinct(batch) |>
    expand_grid(chr = paste0("chr", c(1:22, "X"))) |>
    write_tsv("./array_subset.tsv", col_names = FALSE)

batch_ids |> 
    filter(vcf_id %in% unique(to_fix_df$vcfid)) |>
    write_tsv("./data/metadata.tsv", col_names = FALSE)

batch_ids |> 
    filter(vcf_id %in% unique(to_fix_df$vcfid)) |>
    select(batch, vcf_id) |>
    group_by(batch) |>
    summarise(vcf_id = list(vcf_id)) |>
    ungroup() |>
    mutate(f = paste0("./data/", batch, ".txt")) |>
    {function(dat) walk2(dat$vcf_id, dat$f, ~write_lines(.x, .y))}()

# Fix mapping to FASTQ files
rnaseq_meta <- "/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/B_cells_rnaseq/metadata.tsv" |>
    read_tsv(col_types = c(.default = "c"))

to_fix_df |>
    mutate(bam_sample_id = paste0(donor_bam, ".1"),
	   vcf_sample_id = paste0(donor_vcf, ".1")) |>
    select(subject_id = donor_bam, sample_id = bam_sample_id, stim = bam_stim, 
	   vcf_donor_id = donor_vcf, vcfid, vcf_sample_id) |>
    left_join(rnaseq_meta, join_by(subject_id, sample_id, stim)) |>
    select(subject_id = vcf_donor_id, sample_id = vcf_sample_id, stim, R1, R2) |>
    write_tsv("./data/fix_metadata.tsv", col_names = FALSE)

to_fix_df |>
    mutate(bam_sample_id = paste0(donor_bam, ".1"),
	   vcf_sample_id = paste0(donor_vcf, ".1")) |>
    select(subject_id = donor_bam, sample_id = bam_sample_id, stim = bam_stim, 
	   vcf_donor_id = donor_vcf, vcfid, vcf_sample_id) |>
    left_join(rnaseq_meta, join_by(subject_id, sample_id, stim)) |>
    select(subject_id = vcf_donor_id, sample_id = vcf_sample_id, stim, vcfid) |>
    write_tsv("./array_ase.tsv", col_names = FALSE)
