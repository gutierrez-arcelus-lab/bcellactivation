CITE-seq: Core Processing and Differential Expression
================

This directory contains the scripts for processing the multi-modal
single-cell data, performing rigorous quality control, and conducting
differential expression analysis using pseudobulk methodologies.

**These scripts should be executed in the following order:**

#### Step 1. `run_seurat_v4.R` (R Script)

- **Action:** Processes the multi-modal 10x outputs. It initializes the
  Seurat object, normalizes RNA, ADT, and HTO assays, demultiplexes
  hashing antibodies (`demuxmix`), and integrates the SNP-based doublet
  calls from the previous pipeline step (`demuxlet`). It performs
  Harmony batch correction (controlling for both library run and human
  donor) and generates UMAP projections and cluster markers.
- **Input:** Cell Ranger feature matrices and demuxlet calls.
- **Output:** The final filtered Seurat object (`v4_seurat_qced.rds`),
  coordinates, cluster markers, and QC plots.

#### Step 2. `pseudobulk_dge.R` (R Script)

- **Action:** Performs differential expression testing using **edgeR**.
  It aggregates single-cell counts into biological “pseudobulks” grouped
  by donor and stimulation condition (or cluster). This method prevents
  false positives driven by single-cell independence assumptions and
  allows the statistical model to explicitly control for baseline human
  donor variance (`~ 0 + condition + donor_id`).
- **Input:** The cleaned Seurat object and Cell Ranger gene annotations.
- **Output:** A comprehensive Excel workbook
  (`Supplementary_Data_DGE_citeseq.xlsx`) containing ADT condition
  contrasts, RNA condition contrasts, and RNA cluster markers.

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
    ##  package         * version    date (UTC) lib source
    ##  abind             1.4-5      2016-07-21 [2] CRAN (R 4.1.1)
    ##  backports         1.4.1      2021-12-13 [1] CRAN (R 4.1.2)
    ##  brio              1.1.3      2021-11-30 [2] CRAN (R 4.1.1)
    ##  broom             1.0.0      2022-07-01 [2] CRAN (R 4.1.2)
    ##  cachem            1.0.8      2023-05-01 [1] CRAN (R 4.1.2)
    ##  cellranger        1.1.0      2016-07-27 [2] CRAN (R 4.1.1)
    ##  cli               3.6.4      2025-02-13 [1] CRAN (R 4.1.2)
    ##  cluster           2.1.4      2022-08-22 [2] CRAN (R 4.1.2)
    ##  codetools         0.2-18     2020-11-04 [2] CRAN (R 4.1.1)
    ##  colorspace        2.1-1      2024-07-26 [1] CRAN (R 4.1.2)
    ##  cowplot           1.1.3      2024-01-22 [1] CRAN (R 4.1.2)
    ##  crayon            1.5.2      2022-09-29 [1] CRAN (R 4.1.2)
    ##  data.table        1.15.2     2024-02-29 [1] CRAN (R 4.1.2)
    ##  DBI               1.2.3.9005 2024-07-16 [1] Github (r-dbi/DBI@ba50c6c)
    ##  dbplyr            2.4.0      2023-10-26 [1] CRAN (R 4.1.2)
    ##  deldir            1.0-6      2021-10-23 [2] CRAN (R 4.1.1)
    ##  demuxmix        * 1.5.1      2024-07-24 [1] Github (huklein/demuxmix@4021e00)
    ##  desc              1.4.3      2023-12-10 [1] CRAN (R 4.1.2)
    ##  devtools          2.4.3      2021-11-30 [2] CRAN (R 4.1.1)
    ##  digest            0.6.35     2024-03-11 [1] CRAN (R 4.1.2)
    ##  dplyr           * 1.1.4      2023-11-17 [1] CRAN (R 4.1.2)
    ##  edgeR           * 3.36.0     2021-10-26 [2] Bioconductor
    ##  ellipsis          0.3.2      2021-04-29 [1] CRAN (R 4.1.2)
    ##  evaluate          0.23       2023-11-01 [1] CRAN (R 4.1.2)
    ##  fastmap           1.1.1      2023-02-24 [1] CRAN (R 4.1.2)
    ##  fitdistrplus      1.1-8      2022-03-10 [2] CRAN (R 4.1.1)
    ##  forcats         * 0.5.2      2022-08-19 [2] CRAN (R 4.1.2)
    ##  fs                1.6.3      2023-07-20 [1] CRAN (R 4.1.2)
    ##  future            1.33.1     2023-12-22 [1] CRAN (R 4.1.2)
    ##  future.apply      1.11.1     2023-12-21 [1] CRAN (R 4.1.2)
    ##  gargle            1.5.2      2023-07-20 [1] CRAN (R 4.1.2)
    ##  generics          0.1.3      2022-07-05 [1] CRAN (R 4.1.2)
    ##  ggplot2         * 3.5.1      2024-04-23 [1] CRAN (R 4.1.2)
    ##  ggrepel           0.9.6      2024-09-07 [1] CRAN (R 4.1.2)
    ##  ggridges          0.5.3      2021-01-08 [2] CRAN (R 4.1.1)
    ##  globals           0.16.3     2024-03-08 [1] CRAN (R 4.1.2)
    ##  glue              1.8.0      2024-09-30 [1] CRAN (R 4.1.2)
    ##  goftest           1.2-3      2021-10-07 [2] CRAN (R 4.1.1)
    ##  googledrive       2.0.0      2021-07-08 [2] CRAN (R 4.1.1)
    ##  googlesheets4     1.0.1      2022-08-13 [2] CRAN (R 4.1.2)
    ##  gridExtra         2.3        2017-09-09 [1] CRAN (R 4.1.2)
    ##  gtable            0.3.6      2024-10-25 [1] CRAN (R 4.1.2)
    ##  harmony         * 0.1.1      2022-11-16 [2] Github (immunogenomics/harmony@63ebd73)
    ##  haven             2.5.1      2022-08-22 [2] CRAN (R 4.1.2)
    ##  hms               1.1.3      2023-03-21 [1] CRAN (R 4.1.2)
    ##  htmltools         0.5.7      2023-11-03 [1] CRAN (R 4.1.2)
    ##  htmlwidgets       1.6.4      2023-12-06 [1] CRAN (R 4.1.2)
    ##  httpuv            1.6.14     2024-01-26 [1] CRAN (R 4.1.2)
    ##  httr              1.4.7      2023-08-15 [1] CRAN (R 4.1.2)
    ##  ica               1.0-3      2022-07-08 [2] CRAN (R 4.1.2)
    ##  igraph            2.0.2.9015 2024-03-12 [1] Github (igraph/rigraph@3299d31)
    ##  irlba             2.3.5.1    2022-10-03 [1] CRAN (R 4.1.2)
    ##  jsonlite          1.8.8      2023-12-04 [1] CRAN (R 4.1.2)
    ##  KernSmooth        2.23-20    2021-05-03 [2] CRAN (R 4.1.1)
    ##  knitr             1.45       2023-10-30 [1] CRAN (R 4.1.2)
    ##  later             1.3.2      2023-12-06 [1] CRAN (R 4.1.2)
    ##  lattice           0.20-45    2021-09-22 [2] CRAN (R 4.1.1)
    ##  lazyeval          0.2.2      2019-03-15 [1] CRAN (R 4.1.2)
    ##  leiden            0.4.2      2022-05-09 [2] CRAN (R 4.1.1)
    ##  lifecycle         1.0.4      2023-11-07 [1] CRAN (R 4.1.2)
    ##  limma           * 3.50.3     2022-04-07 [2] Bioconductor
    ##  listenv           0.9.1      2024-01-29 [1] CRAN (R 4.1.2)
    ##  lmtest            0.9-40     2022-03-21 [2] CRAN (R 4.1.1)
    ##  locfit            1.5-9.6    2022-07-11 [2] CRAN (R 4.1.2)
    ##  lubridate         1.8.0      2021-10-07 [2] CRAN (R 4.1.1)
    ##  magrittr          2.0.3      2022-03-30 [1] CRAN (R 4.1.2)
    ##  MASS              7.3-58.1   2022-08-03 [2] CRAN (R 4.1.2)
    ##  Matrix            1.6-3      2023-11-14 [1] CRAN (R 4.1.2)
    ##  matrixStats       1.2.0      2023-12-11 [1] CRAN (R 4.1.2)
    ##  memoise           2.0.1      2021-11-26 [1] CRAN (R 4.1.2)
    ##  mgcv              1.8-40     2022-03-29 [2] CRAN (R 4.1.1)
    ##  mime              0.12       2021-09-28 [1] CRAN (R 4.1.2)
    ##  miniUI            0.1.1.1    2018-05-18 [2] CRAN (R 4.1.1)
    ##  modelr            0.1.9      2022-08-19 [2] CRAN (R 4.1.2)
    ##  munsell           0.5.1      2024-04-01 [1] CRAN (R 4.1.2)
    ##  nlme              3.1-159    2022-08-09 [2] CRAN (R 4.1.2)
    ##  parallelly        1.37.1     2024-02-29 [1] CRAN (R 4.1.2)
    ##  patchwork         1.2.0      2024-01-08 [1] CRAN (R 4.1.2)
    ##  pbapply           1.5-0      2021-09-16 [2] CRAN (R 4.1.1)
    ##  pillar            1.10.2     2025-04-05 [1] CRAN (R 4.1.2)
    ##  pkgbuild          1.4.3      2023-12-10 [1] CRAN (R 4.1.2)
    ##  pkgconfig         2.0.3      2019-09-22 [1] CRAN (R 4.1.2)
    ##  pkgload           1.2.4      2021-11-30 [2] CRAN (R 4.1.1)
    ##  plotly            4.10.4     2024-01-13 [1] CRAN (R 4.1.2)
    ##  plyr              1.8.9      2023-10-02 [1] CRAN (R 4.1.2)
    ##  png               0.1-8      2022-11-29 [1] CRAN (R 4.1.2)
    ##  polyclip          1.10-6     2023-09-27 [1] CRAN (R 4.1.2)
    ##  progressr         0.14.0     2023-08-10 [1] CRAN (R 4.1.2)
    ##  promises          1.2.1      2023-08-10 [1] CRAN (R 4.1.2)
    ##  purrr           * 1.2.2      2026-04-10 [1] CRAN (R 4.1.2)
    ##  R6                2.6.1      2025-02-15 [1] CRAN (R 4.1.2)
    ##  RANN              2.6.1      2019-01-08 [2] CRAN (R 4.1.1)
    ##  RColorBrewer      1.1-3      2022-04-03 [1] CRAN (R 4.1.2)
    ##  Rcpp            * 1.0.14     2025-01-12 [1] CRAN (R 4.1.2)
    ##  RcppAnnoy         0.0.19     2021-07-30 [2] CRAN (R 4.1.1)
    ##  readr           * 2.1.5      2024-01-10 [1] CRAN (R 4.1.2)
    ##  readxl            1.4.1      2022-08-17 [2] CRAN (R 4.1.2)
    ##  remotes           2.4.2      2021-11-30 [2] CRAN (R 4.1.1)
    ##  reprex            2.0.2      2022-08-17 [2] CRAN (R 4.1.2)
    ##  reshape2          1.4.4      2020-04-09 [1] CRAN (R 4.1.2)
    ##  reticulate        1.28       2023-01-27 [2] CRAN (R 4.1.2)
    ##  rgeos             0.5-3      2020-05-08 [2] CRAN (R 4.0.2)
    ##  rlang             1.1.5      2025-01-17 [1] CRAN (R 4.1.2)
    ##  rmarkdown         2.26       2024-03-05 [1] CRAN (R 4.1.2)
    ##  ROCR              1.0-11     2020-05-02 [2] CRAN (R 4.1.1)
    ##  rpart             4.1.16     2022-01-24 [2] CRAN (R 4.1.1)
    ##  rprojroot         2.0.4      2023-11-05 [1] CRAN (R 4.1.2)
    ##  rstudioapi        0.15.0     2023-07-07 [1] CRAN (R 4.1.2)
    ##  Rtsne             0.16       2022-04-17 [2] CRAN (R 4.1.1)
    ##  rvest             1.0.3      2022-08-19 [2] CRAN (R 4.1.2)
    ##  scales            1.3.0      2023-11-28 [1] CRAN (R 4.1.2)
    ##  scattermore       0.8        2022-02-14 [2] CRAN (R 4.1.1)
    ##  sctransform       0.3.4      2022-08-20 [2] CRAN (R 4.1.2)
    ##  sessioninfo       1.2.2      2021-12-06 [2] CRAN (R 4.1.1)
    ##  Seurat          * 4.1.1      2022-05-02 [2] CRAN (R 4.1.1)
    ##  SeuratObject    * 4.1.4      2023-09-26 [1] CRAN (R 4.1.2)
    ##  shiny             1.8.0      2023-11-17 [1] CRAN (R 4.1.2)
    ##  sp                2.1-3      2024-01-30 [1] CRAN (R 4.1.2)
    ##  spatstat.core     2.4-4      2022-05-18 [2] CRAN (R 4.1.2)
    ##  spatstat.data     3.0-0      2022-10-21 [2] CRAN (R 4.1.2)
    ##  spatstat.geom     3.0-3      2022-10-25 [2] CRAN (R 4.1.2)
    ##  spatstat.random   3.0-1      2022-11-03 [2] CRAN (R 4.1.2)
    ##  spatstat.sparse   3.0-0      2022-10-21 [2] CRAN (R 4.1.2)
    ##  spatstat.utils    3.0-5      2024-06-17 [1] CRAN (R 4.1.2)
    ##  stringi           1.8.3      2023-12-11 [1] CRAN (R 4.1.2)
    ##  stringr         * 1.5.1      2023-11-14 [1] CRAN (R 4.1.2)
    ##  survival          3.4-0      2022-08-09 [2] CRAN (R 4.1.2)
    ##  tensor            1.5        2012-05-05 [2] CRAN (R 4.1.1)
    ##  testthat          3.1.4      2022-04-26 [2] CRAN (R 4.1.1)
    ##  tibble          * 3.2.1      2023-03-20 [1] CRAN (R 4.1.2)
    ##  tidyr           * 1.3.1      2024-01-24 [1] CRAN (R 4.1.2)
    ##  tidyselect        1.2.1      2024-03-11 [1] CRAN (R 4.1.2)
    ##  tidyverse       * 1.3.2      2022-07-18 [2] CRAN (R 4.1.2)
    ##  tzdb              0.4.0      2023-05-12 [1] CRAN (R 4.1.2)
    ##  usethis           2.1.6      2022-05-25 [2] CRAN (R 4.1.2)
    ##  uwot              0.1.14     2022-08-22 [2] CRAN (R 4.1.2)
    ##  vctrs             0.6.5      2023-12-01 [1] CRAN (R 4.1.2)
    ##  viridisLite       0.4.2      2023-05-02 [1] CRAN (R 4.1.2)
    ##  withr             3.0.2      2024-10-28 [1] CRAN (R 4.1.2)
    ##  writexl         * 1.5.4      2025-04-15 [1] CRAN (R 4.1.2)
    ##  xfun              0.42       2024-02-08 [1] CRAN (R 4.1.2)
    ##  xml2              1.3.6      2023-12-04 [1] CRAN (R 4.1.2)
    ##  xtable            1.8-4      2019-04-21 [1] CRAN (R 4.1.2)
    ##  yaml              2.3.12     2025-12-10 [1] CRAN (R 4.1.2)
    ##  zoo               1.8-12     2023-04-13 [1] CRAN (R 4.1.2)
    ## 
    ##  [1] /home/ch229163/R/4.1/library
    ##  [2] /programs/biogrids/x86_64-linux/r-pkgs-BIOGRIDS/r4.1-20231122
    ##  [3] /programs/biogrids/x86_64-linux/r/4.1/lib/R/library
    ## 
    ## ──────────────────────────────────────────────────────────────────────────────
