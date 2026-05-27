Transcript Quantification
================

This directory contains the pipeline for building the transcriptome
reference, running the quantification of gene expression, and compiling
the outputs into a finalized count matrix.

#### Step 1. `setup.R`

- **Action:** Generates a sample sheet in the format: three columns
  (sample_id, fastq_1, fastq_2) where the fastq columns are
  comma-separated lists of all fastqs for a sample.
- **Input:** Internal sample sheet file.
- **Output:** `data/metadata.tsv` (Headerless table for Slurm array
  parsing)

#### Step 2. `rsem.slurm` (Slurm Job)

- **Action:** Generates a transcriptome FASTA file from the primary
  assembly.
- **Input:** `GRCh38.primary_assembly.genome.fa` and
  `gencode.v38.primary_assembly.annotation.gtf`
- **Output:** `gencode.v38.primary_assembly.transcripts.fa`

#### Step 3. `salmon_index.slurm` (Slurm Job)

- **Action:** Builds the quasi-mapping index required by Salmon using
  the transcript FASTA generated in Step 1.
- **Input:** `gencode.v38.primary_assembly.transcripts.fa`
- **Output:** Compiled Salmon index.

#### Step 4. `salmon_quant.slurm` (Slurm Array Job)

- **Action:** Quantifies transcript abundance for all samples using
  Salmon.
- **Input:** The sample sheet (`metadata.tsv`) and the Salmon index.
- **Output:** Transcript quantification directories for each sample
  containing `quant.sf` files (Saved in the `/results/` subdirectory).

#### Step 5. `compile_results.R` (R Script)

- **Action:** Imports transcript-level quantifications from the Salmon
  output directories, and aggregates the data into total gene-level read
  counts and TPMs.
- **Input:** `quant.sf` files, `metadata.tsv`, and
  `gencode.v38.primary_assembly.annotation.gtf.gz`.
- **Output:** `expression_pooled_reps.tsv` (Final gene-level counts) and
  `expression_pooled_reps_transcriptlevel.tsv` (Transcript-level
  counts).

#### R Session Information

    ## ─ Session info ───────────────────────────────────────────────────────────────
    ##  setting  value
    ##  version  R version 4.1.2 (2021-11-01)
    ##  os       Rocky Linux 8.10 (Green Obsidian)
    ##  system   x86_64, linux-gnu
    ##  ui       X11
    ##  language (EN)
    ##  collate  en_US.UTF-8
    ##  ctype    en_US.UTF-8
    ##  tz       America/New_York
    ##  date     2026-05-26
    ##  pandoc   2.19.2 @ /programs/biogrids/x86_64-linux/rstudio/2022.02.3/bin// (via rmarkdown)
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────
    ##  package              * version    date (UTC) lib source
    ##  backports              1.4.1      2021-12-13 [1] CRAN (R 4.1.2)
    ##  Biobase                2.54.0     2021-10-26 [1] Bioconductor
    ##  BiocGenerics         * 0.40.0     2021-10-26 [1] Bioconductor
    ##  BiocIO                 1.4.0      2021-10-26 [2] Bioconductor
    ##  BiocParallel           1.28.3     2021-12-09 [2] Bioconductor
    ##  Biostrings             2.62.0     2021-10-26 [1] Bioconductor
    ##  bitops                 1.0-7      2021-04-24 [1] CRAN (R 4.1.2)
    ##  brio                   1.1.3      2021-11-30 [2] CRAN (R 4.1.1)
    ##  broom                  1.0.0      2022-07-01 [2] CRAN (R 4.1.2)
    ##  cachem                 1.0.8      2023-05-01 [1] CRAN (R 4.1.2)
    ##  cellranger             1.1.0      2016-07-27 [2] CRAN (R 4.1.1)
    ##  cli                    3.6.4      2025-02-13 [1] CRAN (R 4.1.2)
    ##  colorspace             2.1-1      2024-07-26 [1] CRAN (R 4.1.2)
    ##  crayon                 1.5.2      2022-09-29 [1] CRAN (R 4.1.2)
    ##  DBI                    1.2.3.9005 2024-07-16 [1] Github (r-dbi/DBI@ba50c6c)
    ##  dbplyr                 2.4.0      2023-10-26 [1] CRAN (R 4.1.2)
    ##  DelayedArray           0.20.0     2021-10-26 [2] Bioconductor
    ##  desc                   1.4.3      2023-12-10 [1] CRAN (R 4.1.2)
    ##  devtools               2.4.3      2021-11-30 [2] CRAN (R 4.1.1)
    ##  digest                 0.6.35     2024-03-11 [1] CRAN (R 4.1.2)
    ##  dplyr                * 1.1.4      2023-11-17 [1] CRAN (R 4.1.2)
    ##  ellipsis               0.3.2      2021-04-29 [1] CRAN (R 4.1.2)
    ##  evaluate               0.23       2023-11-01 [1] CRAN (R 4.1.2)
    ##  fastmap                1.1.1      2023-02-24 [1] CRAN (R 4.1.2)
    ##  forcats              * 0.5.2      2022-08-19 [2] CRAN (R 4.1.2)
    ##  fs                     1.6.3      2023-07-20 [1] CRAN (R 4.1.2)
    ##  gargle                 1.5.2      2023-07-20 [1] CRAN (R 4.1.2)
    ##  generics               0.1.3      2022-07-05 [1] CRAN (R 4.1.2)
    ##  GenomeInfoDb         * 1.30.1     2022-01-30 [2] Bioconductor
    ##  GenomeInfoDbData       1.2.7      2023-07-05 [1] Bioconductor
    ##  GenomicAlignments      1.30.0     2021-10-26 [2] Bioconductor
    ##  GenomicRanges        * 1.46.1     2021-11-18 [2] Bioconductor
    ##  ggplot2              * 3.5.1      2024-04-23 [1] CRAN (R 4.1.2)
    ##  glue                   1.8.0      2024-09-30 [1] CRAN (R 4.1.2)
    ##  googledrive            2.0.0      2021-07-08 [2] CRAN (R 4.1.1)
    ##  googlesheets4          1.0.1      2022-08-13 [2] CRAN (R 4.1.2)
    ##  gtable                 0.3.6      2024-10-25 [1] CRAN (R 4.1.2)
    ##  haven                  2.5.1      2022-08-22 [2] CRAN (R 4.1.2)
    ##  hms                    1.1.3      2023-03-21 [1] CRAN (R 4.1.2)
    ##  htmltools              0.5.7      2023-11-03 [1] CRAN (R 4.1.2)
    ##  httr                   1.4.7      2023-08-15 [1] CRAN (R 4.1.2)
    ##  IRanges              * 2.28.0     2021-10-26 [1] Bioconductor
    ##  jsonlite               1.8.8      2023-12-04 [1] CRAN (R 4.1.2)
    ##  knitr                  1.45       2023-10-30 [1] CRAN (R 4.1.2)
    ##  lattice                0.20-45    2021-09-22 [2] CRAN (R 4.1.1)
    ##  lifecycle              1.0.4      2023-11-07 [1] CRAN (R 4.1.2)
    ##  lubridate              1.8.0      2021-10-07 [2] CRAN (R 4.1.1)
    ##  magrittr               2.0.3      2022-03-30 [1] CRAN (R 4.1.2)
    ##  Matrix                 1.6-3      2023-11-14 [1] CRAN (R 4.1.2)
    ##  MatrixGenerics         1.6.0      2021-10-26 [2] Bioconductor
    ##  matrixStats            1.2.0      2023-12-11 [1] CRAN (R 4.1.2)
    ##  memoise                2.0.1      2021-11-26 [1] CRAN (R 4.1.2)
    ##  modelr                 0.1.9      2022-08-19 [2] CRAN (R 4.1.2)
    ##  munsell                0.5.1      2024-04-01 [1] CRAN (R 4.1.2)
    ##  pillar                 1.10.2     2025-04-05 [1] CRAN (R 4.1.2)
    ##  pkgbuild               1.4.3      2023-12-10 [1] CRAN (R 4.1.2)
    ##  pkgconfig              2.0.3      2019-09-22 [1] CRAN (R 4.1.2)
    ##  pkgload                1.2.4      2021-11-30 [2] CRAN (R 4.1.1)
    ##  purrr                * 1.2.2      2026-04-10 [1] CRAN (R 4.1.2)
    ##  R6                     2.6.1      2025-02-15 [1] CRAN (R 4.1.2)
    ##  RCurl                  1.98-1.14  2024-01-09 [1] CRAN (R 4.1.2)
    ##  readr                * 2.1.5      2024-01-10 [1] CRAN (R 4.1.2)
    ##  readxl                 1.4.1      2022-08-17 [2] CRAN (R 4.1.2)
    ##  remotes                2.4.2      2021-11-30 [2] CRAN (R 4.1.1)
    ##  reprex                 2.0.2      2022-08-17 [2] CRAN (R 4.1.2)
    ##  restfulr               0.0.15     2022-06-16 [2] CRAN (R 4.1.2)
    ##  rjson                  0.2.21     2022-01-09 [2] CRAN (R 4.1.1)
    ##  rlang                  1.1.5      2025-01-17 [1] CRAN (R 4.1.2)
    ##  rmarkdown              2.26       2024-03-05 [1] CRAN (R 4.1.2)
    ##  rprojroot              2.0.4      2023-11-05 [1] CRAN (R 4.1.2)
    ##  Rsamtools              2.10.0     2021-10-26 [2] Bioconductor
    ##  rstudioapi             0.15.0     2023-07-07 [1] CRAN (R 4.1.2)
    ##  rtracklayer          * 1.54.0     2021-10-26 [2] Bioconductor
    ##  rvest                  1.0.3      2022-08-19 [2] CRAN (R 4.1.2)
    ##  S4Vectors            * 0.32.4     2022-03-24 [1] Bioconductor
    ##  scales                 1.3.0      2023-11-28 [1] CRAN (R 4.1.2)
    ##  sessioninfo            1.2.2      2021-12-06 [2] CRAN (R 4.1.1)
    ##  stringi                1.8.3      2023-12-11 [1] CRAN (R 4.1.2)
    ##  stringr              * 1.5.1      2023-11-14 [1] CRAN (R 4.1.2)
    ##  SummarizedExperiment   1.24.0     2021-10-26 [2] Bioconductor
    ##  testthat               3.1.4      2022-04-26 [2] CRAN (R 4.1.1)
    ##  tibble               * 3.2.1      2023-03-20 [1] CRAN (R 4.1.2)
    ##  tidyr                * 1.3.1      2024-01-24 [1] CRAN (R 4.1.2)
    ##  tidyselect             1.2.1      2024-03-11 [1] CRAN (R 4.1.2)
    ##  tidyverse            * 1.3.2      2022-07-18 [2] CRAN (R 4.1.2)
    ##  tximport             * 1.22.0     2021-10-26 [2] Bioconductor
    ##  tzdb                   0.4.0      2023-05-12 [1] CRAN (R 4.1.2)
    ##  usethis                2.1.6      2022-05-25 [2] CRAN (R 4.1.2)
    ##  vctrs                  0.6.5      2023-12-01 [1] CRAN (R 4.1.2)
    ##  withr                  3.0.2      2024-10-28 [1] CRAN (R 4.1.2)
    ##  xfun                   0.42       2024-02-08 [1] CRAN (R 4.1.2)
    ##  XML                    3.99-0.10  2022-06-09 [2] CRAN (R 4.1.2)
    ##  xml2                   1.3.6      2023-12-04 [1] CRAN (R 4.1.2)
    ##  XVector                0.34.0     2021-10-26 [1] Bioconductor
    ##  yaml                   2.3.12     2025-12-10 [1] CRAN (R 4.1.2)
    ##  zlibbioc               1.40.0     2021-10-26 [1] Bioconductor
    ## 
    ##  [1] /home/ch229163/R/4.1/library
    ##  [2] /programs/biogrids/x86_64-linux/r-pkgs-BIOGRIDS/r4.1-20231122
    ##  [3] /programs/biogrids/x86_64-linux/r/4.1/lib/R/library
    ## 
    ## ──────────────────────────────────────────────────────────────────────────────
