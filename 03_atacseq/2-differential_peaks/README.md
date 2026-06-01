ATAC-seq: Differential Peak Accessibility
================

This directory contains the scripts required to evaluate differential
chromatin accessibility between experimental conditions using
**DESeq2**.

#### Step 1. `deseq2.R` (R Script)

- **Action:** Parses the nf-core featureCounts matrix and runs DESeq2.
- **Input:** The `samplesheet.csv` and the consensus peaks count matrix
  (`consensus_peaks.mLb.clN.featureCounts.txt`) generated in the
  `1-processing` directory.
- **Output:**
  - `metadata.tsv` (Cleaned design matrix)
  - `deseq2_data.Rdata` (Raw DESeq2 and regularized log-transformed
    objects)
  - `pcadata_*peaks.rds` (Coordinates for downstream PCA visualizations)
  - `[TREAT]vs[CONTROL].tsv` (Differential accessibility summary tables
    appended with peak genomic coordinates).

------------------------------------------------------------------------

### R Session Information

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
    ##  date     2026-06-01
    ##  pandoc   2.19.2 @ /programs/biogrids/x86_64-linux/rstudio/2022.02.3/bin// (via rmarkdown)
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────
    ##  package              * version    date (UTC) lib source
    ##  affy                   1.72.0     2021-10-26 [2] Bioconductor
    ##  affyio                 1.64.0     2021-10-26 [2] Bioconductor
    ##  annotate               1.72.0     2021-10-26 [2] Bioconductor
    ##  AnnotationDbi          1.56.2     2021-11-09 [1] Bioconductor
    ##  backports              1.4.1      2021-12-13 [1] CRAN (R 4.1.2)
    ##  Biobase              * 2.54.0     2021-10-26 [1] Bioconductor
    ##  BiocGenerics         * 0.40.0     2021-10-26 [1] Bioconductor
    ##  BiocManager            1.30.25    2024-08-28 [1] CRAN (R 4.1.2)
    ##  BiocParallel         * 1.28.3     2021-12-09 [2] Bioconductor
    ##  Biostrings             2.62.0     2021-10-26 [1] Bioconductor
    ##  bit                    4.0.5      2022-11-15 [1] CRAN (R 4.1.2)
    ##  bit64                  4.0.5      2020-08-30 [1] CRAN (R 4.1.2)
    ##  bitops                 1.0-7      2021-04-24 [1] CRAN (R 4.1.2)
    ##  blob                   1.2.4      2023-03-17 [1] CRAN (R 4.1.2)
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
    ##  DESeq2               * 1.34.0     2021-10-26 [2] Bioconductor
    ##  devtools               2.4.3      2021-11-30 [2] CRAN (R 4.1.1)
    ##  digest                 0.6.35     2024-03-11 [1] CRAN (R 4.1.2)
    ##  dplyr                * 1.1.4      2023-11-17 [1] CRAN (R 4.1.2)
    ##  ellipsis               0.3.2      2021-04-29 [1] CRAN (R 4.1.2)
    ##  evaluate               0.23       2023-11-01 [1] CRAN (R 4.1.2)
    ##  fastmap                1.1.1      2023-02-24 [1] CRAN (R 4.1.2)
    ##  forcats              * 0.5.2      2022-08-19 [2] CRAN (R 4.1.2)
    ##  fs                     1.6.3      2023-07-20 [1] CRAN (R 4.1.2)
    ##  gargle                 1.5.2      2023-07-20 [1] CRAN (R 4.1.2)
    ##  genefilter             1.76.0     2021-10-26 [2] Bioconductor
    ##  geneplotter            1.72.0     2021-10-26 [2] Bioconductor
    ##  generics               0.1.3      2022-07-05 [1] CRAN (R 4.1.2)
    ##  GenomeInfoDb         * 1.30.1     2022-01-30 [2] Bioconductor
    ##  GenomeInfoDbData       1.2.7      2023-07-05 [1] Bioconductor
    ##  GenomicRanges        * 1.46.1     2021-11-18 [2] Bioconductor
    ##  ggplot2              * 3.5.1      2024-04-23 [1] CRAN (R 4.1.2)
    ##  glue                 * 1.8.0      2024-09-30 [1] CRAN (R 4.1.2)
    ##  googledrive            2.0.0      2021-07-08 [2] CRAN (R 4.1.1)
    ##  googlesheets4          1.0.1      2022-08-13 [2] CRAN (R 4.1.2)
    ##  gtable                 0.3.6      2024-10-25 [1] CRAN (R 4.1.2)
    ##  haven                  2.5.1      2022-08-22 [2] CRAN (R 4.1.2)
    ##  hms                    1.1.3      2023-03-21 [1] CRAN (R 4.1.2)
    ##  htmltools              0.5.7      2023-11-03 [1] CRAN (R 4.1.2)
    ##  httr                   1.4.7      2023-08-15 [1] CRAN (R 4.1.2)
    ##  IRanges              * 2.28.0     2021-10-26 [1] Bioconductor
    ##  jsonlite               1.8.8      2023-12-04 [1] CRAN (R 4.1.2)
    ##  KEGGREST               1.34.0     2021-10-26 [1] Bioconductor
    ##  knitr                  1.45       2023-10-30 [1] CRAN (R 4.1.2)
    ##  lattice                0.20-45    2021-09-22 [2] CRAN (R 4.1.1)
    ##  lifecycle              1.0.4      2023-11-07 [1] CRAN (R 4.1.2)
    ##  limma                  3.50.3     2022-04-07 [2] Bioconductor
    ##  locfit                 1.5-9.6    2022-07-11 [2] CRAN (R 4.1.2)
    ##  lubridate              1.8.0      2021-10-07 [2] CRAN (R 4.1.1)
    ##  magrittr               2.0.3      2022-03-30 [1] CRAN (R 4.1.2)
    ##  Matrix                 1.6-3      2023-11-14 [1] CRAN (R 4.1.2)
    ##  MatrixGenerics       * 1.6.0      2021-10-26 [2] Bioconductor
    ##  matrixStats          * 1.2.0      2023-12-11 [1] CRAN (R 4.1.2)
    ##  memoise                2.0.1      2021-11-26 [1] CRAN (R 4.1.2)
    ##  modelr                 0.1.9      2022-08-19 [2] CRAN (R 4.1.2)
    ##  munsell                0.5.1      2024-04-01 [1] CRAN (R 4.1.2)
    ##  pillar                 1.10.2     2025-04-05 [1] CRAN (R 4.1.2)
    ##  pkgbuild               1.4.3      2023-12-10 [1] CRAN (R 4.1.2)
    ##  pkgconfig              2.0.3      2019-09-22 [1] CRAN (R 4.1.2)
    ##  pkgload                1.2.4      2021-11-30 [2] CRAN (R 4.1.1)
    ##  png                    0.1-8      2022-11-29 [1] CRAN (R 4.1.2)
    ##  preprocessCore         1.56.0     2021-10-26 [1] Bioconductor
    ##  purrr                * 1.2.2      2026-04-10 [1] CRAN (R 4.1.2)
    ##  R6                     2.6.1      2025-02-15 [1] CRAN (R 4.1.2)
    ##  RColorBrewer           1.1-3      2022-04-03 [1] CRAN (R 4.1.2)
    ##  Rcpp                   1.0.14     2025-01-12 [1] CRAN (R 4.1.2)
    ##  RCurl                  1.98-1.14  2024-01-09 [1] CRAN (R 4.1.2)
    ##  readr                * 2.1.5      2024-01-10 [1] CRAN (R 4.1.2)
    ##  readxl                 1.4.1      2022-08-17 [2] CRAN (R 4.1.2)
    ##  remotes                2.4.2      2021-11-30 [2] CRAN (R 4.1.1)
    ##  reprex                 2.0.2      2022-08-17 [2] CRAN (R 4.1.2)
    ##  rlang                  1.1.5      2025-01-17 [1] CRAN (R 4.1.2)
    ##  rmarkdown              2.26       2024-03-05 [1] CRAN (R 4.1.2)
    ##  rprojroot              2.0.4      2023-11-05 [1] CRAN (R 4.1.2)
    ##  RSQLite                2.3.5      2024-01-21 [1] CRAN (R 4.1.2)
    ##  rstudioapi             0.15.0     2023-07-07 [1] CRAN (R 4.1.2)
    ##  rvest                  1.0.3      2022-08-19 [2] CRAN (R 4.1.2)
    ##  S4Vectors            * 0.32.4     2022-03-24 [1] Bioconductor
    ##  scales                 1.3.0      2023-11-28 [1] CRAN (R 4.1.2)
    ##  sessioninfo            1.2.2      2021-12-06 [2] CRAN (R 4.1.1)
    ##  stringi                1.8.3      2023-12-11 [1] CRAN (R 4.1.2)
    ##  stringr              * 1.5.1      2023-11-14 [1] CRAN (R 4.1.2)
    ##  SummarizedExperiment * 1.24.0     2021-10-26 [2] Bioconductor
    ##  survival               3.4-0      2022-08-09 [2] CRAN (R 4.1.2)
    ##  testthat               3.1.4      2022-04-26 [2] CRAN (R 4.1.1)
    ##  tibble               * 3.2.1      2023-03-20 [1] CRAN (R 4.1.2)
    ##  tidyr                * 1.3.1      2024-01-24 [1] CRAN (R 4.1.2)
    ##  tidyselect             1.2.1      2024-03-11 [1] CRAN (R 4.1.2)
    ##  tidyverse            * 1.3.2      2022-07-18 [2] CRAN (R 4.1.2)
    ##  tzdb                   0.4.0      2023-05-12 [1] CRAN (R 4.1.2)
    ##  usethis                2.1.6      2022-05-25 [2] CRAN (R 4.1.2)
    ##  vctrs                  0.6.5      2023-12-01 [1] CRAN (R 4.1.2)
    ##  vsn                  * 3.62.0     2021-10-26 [2] Bioconductor
    ##  withr                  3.0.2      2024-10-28 [1] CRAN (R 4.1.2)
    ##  xfun                   0.42       2024-02-08 [1] CRAN (R 4.1.2)
    ##  XML                    3.99-0.10  2022-06-09 [2] CRAN (R 4.1.2)
    ##  xml2                   1.3.6      2023-12-04 [1] CRAN (R 4.1.2)
    ##  xtable                 1.8-4      2019-04-21 [1] CRAN (R 4.1.2)
    ##  XVector                0.34.0     2021-10-26 [1] Bioconductor
    ##  yaml                   2.3.12     2025-12-10 [1] CRAN (R 4.1.2)
    ##  zlibbioc               1.40.0     2021-10-26 [1] Bioconductor
    ## 
    ##  [1] /home/ch229163/R/4.1/library
    ##  [2] /programs/biogrids/x86_64-linux/r-pkgs-BIOGRIDS/r4.1-20231122
    ##  [3] /programs/biogrids/x86_64-linux/r/4.1/lib/R/library
    ## 
    ## ──────────────────────────────────────────────────────────────────────────────
