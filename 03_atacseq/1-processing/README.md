### ATAC-seq Processing 

This directory contains the setup and execution scripts for the primary ATAC-seq data processing pipeline using **[nf-core/atacseq v2.0](https://nf-co.re/atacseq/2.0)**. The pipeline handles adapter trimming, alignment to the human genome, mitochondrial read filtering, duplicate removal, and narrow peak calling via MACS2.

#### Prerequisite: Sample Sheet (`samplesheet.csv`)
The pipeline requires a standard nf-core formatted CSV file. It must contain headers for `sample`, `fastq_1`, `fastq_2`, and `replicate`. 

#### Step 1. `launch.slurm` (Slurm Job)
* **Action:** Acts as the orchestrator node for Nextflow. It sets up the required temporary caching directories for Singularity containers and executes the pipeline using our custom cluster configuration. 
* **Input:** `samplesheet.csv`, `e2.config` (custom executor settings), and the ENCODE hg38 v3 blacklist BED file.
* **Output:** Processed alignments (BAMs), called peaks, and comprehensive MultiQC quality control reports located in the `./results/` directory.

**Note for External Users:** The `e2.config` file contains institutional-specific Slurm executor limits (e.g., queues, max memory, max time). If running this outside our cluster, you should provide your own `-profile` (like `docker`, `singularity`, or a custom institutional profile) instead of `e2.config`.
