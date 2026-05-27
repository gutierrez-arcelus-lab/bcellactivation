# ==============================================================================
# Purpose: Parse RNA-seq fastq filenames to generate a structured metadata
#          sample sheet (`metadata.tsv`). This sample sheet is used as the 
#          input array for downstream SLURM processing scripts (e.g., FastQC).
#
# NOTE FOR EXTERNAL USERS: 
# The file paths hardcoded in `fastq_files` below reflect the internal server 
# directory structure of the Gutierrez-Arcelus lab. If you are reproducing this 
# analysis using data downloaded from dbGaP, please point the `fastq_files` 
# object to your local directory containing the downloaded .fastq.gz files.
# For example: list.files("./dbgap_data/", full.names = TRUE, pattern = "\\.fastq\\.gz")
# ==============================================================================

library(tidyverse)

# 1. Locate all raw FASTQ files
fastq_files <-
    file.path("./data/B_cells_rnaseq", 
	      c("34.198.31.178/221025_MG10430_fastq", 
		"resequencing/231227_MG10430_fastq",
		"resequencing/240122_MG10430_fastq",
		"resequencing/240320_MG10430_fastq")) |>
    map(~list.files(., full.names = TRUE, pattern = "\\.fastq\\.gz")) |>
    unlist()

# 2. Parse filenames into a structured metadata table
meta <- 
    tibble(f = fastq_files) |>
    # Extract the base filename and the read direction (R1 or R2)
    mutate(basen = basename(f),
	   r = str_extract(basen, "_(R[12])_", group = 1)) |>
    # Extract metadata from the complex naming convention using Regex
    # Pattern matches: {batch}_{dummy}_{donor_id}_{replic_id}_{stim}
    extract(basen, 
	    c("batch", "dummy", "donor_id", "replic_id", "stim"),
	    "(\\d+)_(MG_)?(\\d+)[_-](\\d[_-])?([^_-]+)",
	    remove = FALSE) |>
    select(batch, donor_id, replic_id, stim, f, basen, r) |>
    # Clean up formatting
    mutate(replic_id = parse_number(replic_id),
	   replic_id = replace_na(replic_id, 1),
	   stim = recode(stim, "TR7" = "TLR7")) |>
    # Ensure we only keep samples that have both an R1 and R2 file
    group_by(batch, donor_id, replic_id, stim) |>
    filter(all(c("R1", "R2") %in% r)) |>
    ungroup() |>
    # Handle libraries with resequencing to increase depth: 
    # collapse their file paths into a comma-separated string for easy SLURM passing
    summarise(fastq = paste(f, collapse = ","),
	      .by = c(donor_id, replic_id, stim, r)) |>
    # Pivot to wide format so each sample has one row with an 'R1' and 'R2' column
    pivot_wider(names_from = r, values_from = fastq)

# 3. Export the final sample sheet for the SLURM arrays
write_tsv(meta, "./data/metadata.tsv")
