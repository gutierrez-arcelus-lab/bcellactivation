### Transcript Quantification 

This directory contains the Slurm scripts used to build a transcriptome reference, generate a decoy-aware Salmon index, and quantify transcript-level expression from the adapter-trimmed RNA-seq reads.

#### Prerequisite: Sample Sheet (`metadata.tsv`)
The quantification step utilizes the same `metadata.tsv` sample sheet generated/formatted during the mapping step. Ensure this file is present in the `data/` subdirectory.

---

**These scripts should be executed in the following order:**

#### Step 1. `rsem.slurm` (Slurm Job)
* **Action:** Parses the GRCh38 primary assembly genome and GENCODE v41 annotations to extract a transcriptome FASTA file. 
* **Input:** `GRCh38.primary_assembly.genome.fa` and `gencode.v41.primary_assembly.annotation.gtf`.
* **Output:** `gencode.v41.primary_assembly.transcripts.fa` (Used as input for Salmon indexing).

#### Step 2. `salmon_index.slurm` (Slurm Job)
* **Action:** Builds a decoy-aware quasi-mapping index for Salmon. 
* **Input:** The genome FASTA and the `transcripts.fa` generated in Step 1.
* **Output:** Compiled decoy-aware Salmon index inside the `./index/` directory.

#### Step 3. `salmon_quant.slurm` (Slurm Array Job)
* **Action:** Quantifies transcript abundance for each sample using Salmon. 
* **Input:** The decoy-aware Salmon index and the trimmed paired-end `.fq.gz` files (paths provided via `metadata.tsv`).
* **Output:** Transcript-level quantification files (`quant.sf`) stored in the `./results/` directory. 
