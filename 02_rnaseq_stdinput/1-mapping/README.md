### Read Mapping & Junction Extraction

This directory contains the Slurm scripts used to align the adapter-trimmed RNA-seq reads to the human reference genome and extract splice junctions for downstream splicing analyses.

#### Prerequisite: Mapping Sample Sheet (`metadata.tsv`)
Before running the mapping array, ensure you have a `metadata.tsv` file in the `data/` subdirectory. Unlike the long-format sheet used for QC, if a sample has multiple trimmed FASTQ files, they should be collapsed into a single row with their paths separated by commas (no spaces).

**Expected Format (3 columns, with header):**
| sample_id | fastq_R1 | fastq_R2 |
| :--- | :--- | :--- |
| 10430_1_CD40L | `/path/to/S1_L001_R1_val_1.fq.gz,/path/to/S1_L002_R1_val_1.fq.gz` | `/path/to/S1_L001_R2_val_2.fq.gz,/path/to/S1_L002_R2_val_2.fq.gz` |
| 10431_1_TLR9 | `/path/to/S2_L001_R1_val_1.fq.gz` | `/path/to/S2_L001_R2_val_2.fq.gz` |

---

**These scripts should be executed in the following order:**

#### Step 1. `star_index.slurm` (Slurm Job)
* **Action:** Generates the STAR genome index. 
* **Input:** `GRCh38.primary_assembly.genome.fa` and `gencode.v41.primary_assembly.annotation.gtf`.
* **Output:** Compiled STAR genome index saved in the `./index` directory.

#### Step 2. `star_map.slurm` (Slurm Array Job)
* **Action:** Aligns the adapter-trimmed reads to the genome using STAR's 2-pass mode. Following alignment, it isolates strictly uniquely mapping reads (MAPQ 255) on the primary chromosomes (1-22, X) and uses `regtools` to extract splice junctions for downstream splicing analysis tools.
* **Input:** The `metadata.tsv` sample sheet, the compiled STAR index, and the trimmed `.fq.gz` files.
* **Output:** Coordinate-sorted BAM files (`*Aligned.sortedByCoord.out.bam`) and extracted splice junction files (`*.junc`)
