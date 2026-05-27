Differential Splicing Analysis
================

This directory contains the pipeline for identifying alternative
splicing events using **LeafCutter**. The workflow takes the
high-confidence splice junctions extracted during the STAR mapping step,
clusters overlapping introns, and performs differential excision
analysis across stimulation conditions.

**These scripts should be executed in the following order:**

#### Step 1. `setup.R` (R Script)

- **Action:** Parses the main `metadata.tsv` to dynamically generate the
  pairwise contrast files required by LeafCutter.
- **Input:** `../1-mapping/data/metadata.tsv`
- **Output:** `junction_files.txt`, `contrasts.txt`, and sample grouping
  `.tsv` files for each pairwise comparison.

#### Step 2. `cluster_introns.sh` (Bash Script)

- **Action:** Runs LeafCutter’s core Python clustering algorithm. It
  aggregates overlapping split reads (junctions) across all samples into
  distinct clusters representing alternative splicing events, discarding
  rare junctions.
- **Input:** `junction_files.txt` and the raw `.junc` files from the
  mapping directory.
- **Output:** `clusters_perind_numers.counts.gz` (Cluster count matrix).
- *Note: This step requires a Python 2.7 environment and a local
  installation of LeafCutter.*

#### Step 3. `make_exon_file.sh` (Bash Script)

- **Action:** Parses the GENCODE annotation file to generate a
  compressed exon definitions file. This allows LeafCutter to map splice
  junctions back to known genes and exons.
- **Input:** `gencode.v41.primary_assembly.annotation.gtf`
- **Output:** `exon.txt.gz`

#### Step 4. `make_leafviz_annot.sh` (Bash Script)

- **Action:** Builds the custom genomic and transcriptomic annotation
  database required by the LeafViz application.
- **Input:** `gencode.v41.primary_assembly.annotation.gtf`
- **Output:** Annotation database generated in the
  `./data/leafviz_annot/` directory.

#### Step 5. `leafcutter_ds.slurm` (Slurm Array Job)

- **Action:** Uses a Dirichlet-multinomial generalized linear model to
  test for differential intron excision between biological conditions.
  The Slurm array iterates through all pairwise contrasts defined in
  Step 1.
- **Input:** The count matrix, contrast group mappings, and exon
  definitions.
- **Output:** Statistical results (`_cluster_significance.txt`) and
  effect sizes (`_effect_sizes.txt`) stored in condition-specific
  folders within `./results/`.

#### Step 6. `prepare_leafviz.slurm` (Slurm Array Job)

- **Action:** Compiles the statistical outputs into `.Rdata` objects
  required for the interactive LeafViz visualization application.
- **Input:** The differential splicing outputs and pre-built LeafViz
  annotation files.
- **Output:** `[CONTRAST].Rdata` visualization objects.

#### Step 7. `run_leafviz.R` (Local Utility Script)

- **Action:** Launches the LeafViz interactive Shiny server locally,
  allowing the user to browse significant splicing clusters, view
  differential splicing boxplots, and explore gene tracks.
- **Input:** The compiled `.Rdata` files from Step 6.

------------------------------------------------------------------------

### Session Information

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
    ##  date     2026-05-27
    ##  pandoc   2.19.2 @ /programs/biogrids/x86_64-linux/rstudio/2022.02.3/bin// (via rmarkdown)
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────
    ##  package       * version    date (UTC) lib source
    ##  backports       1.4.1      2021-12-13 [1] CRAN (R 4.1.2)
    ##  brio            1.1.3      2021-11-30 [2] CRAN (R 4.1.1)
    ##  broom           1.0.0      2022-07-01 [2] CRAN (R 4.1.2)
    ##  cachem          1.0.8      2023-05-01 [1] CRAN (R 4.1.2)
    ##  cellranger      1.1.0      2016-07-27 [2] CRAN (R 4.1.1)
    ##  cli             3.6.4      2025-02-13 [1] CRAN (R 4.1.2)
    ##  colorspace      2.1-1      2024-07-26 [1] CRAN (R 4.1.2)
    ##  crayon          1.5.2      2022-09-29 [1] CRAN (R 4.1.2)
    ##  DBI             1.2.3.9005 2024-07-16 [1] Github (r-dbi/DBI@ba50c6c)
    ##  dbplyr          2.4.0      2023-10-26 [1] CRAN (R 4.1.2)
    ##  desc            1.4.3      2023-12-10 [1] CRAN (R 4.1.2)
    ##  devtools        2.4.3      2021-11-30 [2] CRAN (R 4.1.1)
    ##  digest          0.6.35     2024-03-11 [1] CRAN (R 4.1.2)
    ##  dplyr         * 1.1.4      2023-11-17 [1] CRAN (R 4.1.2)
    ##  ellipsis        0.3.2      2021-04-29 [1] CRAN (R 4.1.2)
    ##  evaluate        0.23       2023-11-01 [1] CRAN (R 4.1.2)
    ##  fastmap         1.1.1      2023-02-24 [1] CRAN (R 4.1.2)
    ##  forcats       * 0.5.2      2022-08-19 [2] CRAN (R 4.1.2)
    ##  fs              1.6.3      2023-07-20 [1] CRAN (R 4.1.2)
    ##  gargle          1.5.2      2023-07-20 [1] CRAN (R 4.1.2)
    ##  generics        0.1.3      2022-07-05 [1] CRAN (R 4.1.2)
    ##  ggplot2       * 3.5.1      2024-04-23 [1] CRAN (R 4.1.2)
    ##  ggrepel         0.9.6      2024-09-07 [1] CRAN (R 4.1.2)
    ##  glue          * 1.8.0      2024-09-30 [1] CRAN (R 4.1.2)
    ##  googledrive     2.0.0      2021-07-08 [2] CRAN (R 4.1.1)
    ##  googlesheets4   1.0.1      2022-08-13 [2] CRAN (R 4.1.2)
    ##  gridExtra       2.3        2017-09-09 [1] CRAN (R 4.1.2)
    ##  gtable          0.3.6      2024-10-25 [1] CRAN (R 4.1.2)
    ##  haven           2.5.1      2022-08-22 [2] CRAN (R 4.1.2)
    ##  hms             1.1.3      2023-03-21 [1] CRAN (R 4.1.2)
    ##  htmltools       0.5.7      2023-11-03 [1] CRAN (R 4.1.2)
    ##  httr            1.4.7      2023-08-15 [1] CRAN (R 4.1.2)
    ##  jsonlite        1.8.8      2023-12-04 [1] CRAN (R 4.1.2)
    ##  knitr           1.45       2023-10-30 [1] CRAN (R 4.1.2)
    ##  leafviz       * 0.1.0      2025-05-27 [1] Github (jackhump/leafviz@0ba9688)
    ##  lifecycle       1.0.4      2023-11-07 [1] CRAN (R 4.1.2)
    ##  lubridate       1.8.0      2021-10-07 [2] CRAN (R 4.1.1)
    ##  magrittr        2.0.3      2022-03-30 [1] CRAN (R 4.1.2)
    ##  memoise         2.0.1      2021-11-26 [1] CRAN (R 4.1.2)
    ##  modelr          0.1.9      2022-08-19 [2] CRAN (R 4.1.2)
    ##  munsell         0.5.1      2024-04-01 [1] CRAN (R 4.1.2)
    ##  pillar          1.10.2     2025-04-05 [1] CRAN (R 4.1.2)
    ##  pkgbuild        1.4.3      2023-12-10 [1] CRAN (R 4.1.2)
    ##  pkgconfig       2.0.3      2019-09-22 [1] CRAN (R 4.1.2)
    ##  pkgload         1.2.4      2021-11-30 [2] CRAN (R 4.1.1)
    ##  purrr         * 1.2.2      2026-04-10 [1] CRAN (R 4.1.2)
    ##  R6              2.6.1      2025-02-15 [1] CRAN (R 4.1.2)
    ##  Rcpp            1.0.14     2025-01-12 [1] CRAN (R 4.1.2)
    ##  readr         * 2.1.5      2024-01-10 [1] CRAN (R 4.1.2)
    ##  readxl          1.4.1      2022-08-17 [2] CRAN (R 4.1.2)
    ##  remotes         2.4.2      2021-11-30 [2] CRAN (R 4.1.1)
    ##  reprex          2.0.2      2022-08-17 [2] CRAN (R 4.1.2)
    ##  rlang           1.1.5      2025-01-17 [1] CRAN (R 4.1.2)
    ##  rmarkdown       2.26       2024-03-05 [1] CRAN (R 4.1.2)
    ##  rprojroot       2.0.4      2023-11-05 [1] CRAN (R 4.1.2)
    ##  rstudioapi      0.15.0     2023-07-07 [1] CRAN (R 4.1.2)
    ##  rvest           1.0.3      2022-08-19 [2] CRAN (R 4.1.2)
    ##  scales          1.3.0      2023-11-28 [1] CRAN (R 4.1.2)
    ##  sessioninfo     1.2.2      2021-12-06 [2] CRAN (R 4.1.1)
    ##  stringi         1.8.3      2023-12-11 [1] CRAN (R 4.1.2)
    ##  stringr       * 1.5.1      2023-11-14 [1] CRAN (R 4.1.2)
    ##  testthat        3.1.4      2022-04-26 [2] CRAN (R 4.1.1)
    ##  tibble        * 3.2.1      2023-03-20 [1] CRAN (R 4.1.2)
    ##  tidyr         * 1.3.1      2024-01-24 [1] CRAN (R 4.1.2)
    ##  tidyselect      1.2.1      2024-03-11 [1] CRAN (R 4.1.2)
    ##  tidyverse     * 1.3.2      2022-07-18 [2] CRAN (R 4.1.2)
    ##  tzdb            0.4.0      2023-05-12 [1] CRAN (R 4.1.2)
    ##  usethis         2.1.6      2022-05-25 [2] CRAN (R 4.1.2)
    ##  vctrs           0.6.5      2023-12-01 [1] CRAN (R 4.1.2)
    ##  withr           3.0.2      2024-10-28 [1] CRAN (R 4.1.2)
    ##  xfun            0.42       2024-02-08 [1] CRAN (R 4.1.2)
    ##  xml2            1.3.6      2023-12-04 [1] CRAN (R 4.1.2)
    ##  yaml            2.3.12     2025-12-10 [1] CRAN (R 4.1.2)
    ## 
    ##  [1] /home/ch229163/R/4.1/library
    ##  [2] /programs/biogrids/x86_64-linux/r-pkgs-BIOGRIDS/r4.1-20231122
    ##  [3] /programs/biogrids/x86_64-linux/r/4.1/lib/R/library
    ## 
    ## ──────────────────────────────────────────────────────────────────────────────
