# ==============================================================================
# Description:  Parses MultiQC/FastQC outputs to identify and filter RNA-seq 
#               libraries. Sums reads across lanes/replicates and retains only 
#               samples with >= 2 million total reads.
# Input:        description_hsapiens_qc.tsv (Raw internal metadata)
#               mqc_fastqc_sequence_counts_plot_1.txt (MultiQC output)
# Output:       samples_pass.tsv (List of samples passing the 2M read threshold)
#
# NOTE FOR EXTERNAL REPRODUCIBILITY:
# dbGaP only contains the samples with successfull libaries, therefore external
# users should not run this script.
# ==============================================================================

library(tidyverse)

# Fetch the system environment variable for the working directory
# This is where the output of FASTQ was stored
temp_dir <- Sys.getenv("TEMP_WORK")

# ------------------------------------------------------------------------------
# 1. Parse Metadata
# ------------------------------------------------------------------------------
meta_long <- 
    "../data/description_hsapiens_qc.tsv" |>
    read_tsv() |>
    # Reshape so Read 1 and Read 2 are in the same column
    pivot_longer(fq1:fq2, names_to = "dummy", values_to = "fastq") |>
    # Expand comma-separated FASTQ lists (from pooled lanes) into individual rows
    separate_rows(fastq, sep = ",") |>
    mutate(
	   # parse fastq names
	   fastq = basename(fastq),
           fastq = sub("\\.fq\\.gz", "", fastq),
	   # Extract the lane number using regex
	   lane = sub("^(\\d)_.+$", "\\1", fastq)) |>
    select(-barcode_seq, -dummy)

# ------------------------------------------------------------------------------
# 2. Load MultiQC Sequence Counts
# ------------------------------------------------------------------------------
fastqc <-
    file.path(temp_dir, "fastq/lowinput/multiqc_data/mqc_fastqc_sequence_counts_plot_1.txt") |>
    read_tsv()

# ------------------------------------------------------------------------------
# 3. Merge and Format Data
# ------------------------------------------------------------------------------
fastqc_df <- 
    # Join MultiQC data to the metadata based on the FASTQ filename
    left_join(meta_long, fastqc, join_by(fastq == Sample)) |>
    pivot_longer(`Unique Reads`:`Duplicate Reads`, names_to = "read_type", values_to = "n") |>
    mutate(
	   # Clean up column names and identify if the file is Read 1 or Read 2
	   read_type = sub("\\sReads$", "", read_type),
           fastq = sub("^.+([12])$", "\\1", fastq),
	   # Create an ID containing sample, stim, time, lane, and read number
           id = paste(sample_id, stim, time, paste0("L", lane), paste0("fq", fastq), sep = "_"))

# ------------------------------------------------------------------------------
# 4. Apply 2 Million Read Threshold
# ------------------------------------------------------------------------------
keep_samples <- 
    fastqc_df |>
    # Filter for Read 1 only to prevent double-counting paired-end reads, 
    # and exclude the blank/control samples
    filter(fastq == 1, sample_id != "BLANK") |>
    # Sum the reads across all lanes/technical replicates for each biological sample
    group_by(sample_id, stim, time) |>
    mutate(n = sum(n)) |>
    ungroup() |>
    # Keep only samples with at least 2 million reads
    filter(n >= 2e6) |>
    # Collapse back to unique samples
    distinct(plate, well, sample_id, stim, time) |>
    # Create a plate_well ID to map back to the physical source
    unite("id", c(plate, well), sep = "_", remove = FALSE)

# ------------------------------------------------------------------------------
# 5. Export Final Passed Samples
# ------------------------------------------------------------------------------
keep_samples |>
    # Clean up technical replicate suffixes from the final IDs
    mutate(sample_id = sub("\\.rep\\.\\d", "", sample_id)) |>
    unite("sample_name", c("sample_id", "stim", "time"), sep = "_") |>
    select(sample_id = id, sample_name) |>
    write_tsv("./data/samples_pass.tsv")
