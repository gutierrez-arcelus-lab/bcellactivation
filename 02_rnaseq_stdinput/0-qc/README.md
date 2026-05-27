### Quality Control

This directory contains the Slurm scripts used to prepare standard bulk RNA-seq samples for analysis. The workflow removes adapter sequences, filters for minimum read length, and performs initial quality assessment.

#### Prerequisite: Sample Sheet (`metadata.tsv`)
Before running the QC array, users downloading raw data from dbGaP must create a tab-separated file named `metadata.tsv` inside a `data/` subdirectory. This file dictates the Slurm array execution. 

The file must be in **long format**. If a single biological sample was sequenced across multiple lanes or batches to increase depth, each FASTQ pair must have its own distinct row. 

**Expected Format (3 columns, with header):**
| sample_id | fastq_R1 | fastq_R2 |
| :--- | :--- | :--- |
| 10430_1_CD40L | `/path/to/S1_L001_R1.fq.gz` | `/path/to/S1_L001_R2.fq.gz` |
| 10430_1_CD40L | `/path/to/S1_L002_R1.fq.gz` | `/path/to/S1_L002_R2.fq.gz` |
| 10431_1_TLR9 | `/path/to/S2_L001_R1.fq.gz` | `/path/to/S2_L001_R2.fq.gz` |

#### Step 1. `run_qc.slurm` (Slurm Array Job)
* **Action:** For each FASTQ pair, it runs `trim_galore` to remove adapter sequences and discard reads that drop below 75bp. It automatically runs `fastqc` on the resulting trimmed files.
* **Input:** `metadata.tsv` and raw paired-end `.fastq.gz` files.
* **Output:** Trimmed `.fq.gz` files and `.html`/`.zip` FastQC reports generated in the specified output directory.
