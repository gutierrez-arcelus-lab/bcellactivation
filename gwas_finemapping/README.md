GWAS Fine-mapping Pipeline (TNFSF4 Locus)
================

This directory contains the pipeline for performing statistical
fine-mapping of the systemic lupus erythematosus (SLE) GWAS locus at
*TNFSF4* using SuSiE (Sum of Single Effects).

## Order of Execution

1.  **`setup.R`**: Parses the raw SLE GWAS supplementary tables from the
    literature (Langefeld et al.), isolates the *TNFSF4* locus genomic
    window, and formats the summary statistics for downstream
    processing.
2.  **`munge_lift.sh`**: A shell script that processes the GWAS summary
    statistics using `bcftools +munge`, normalizes them against the
    GRCh37 reference, and generates the liftover versions necessary for
    downstream plotting.
3.  **`infer_ld.sh`**: A shell script that extracts the corresponding
    *TNFSF4* region from the 1000 Genomes Phase 3 European reference
    panel, normalizes the variants, and uses `plink2` to generate the
    correlation matrix.
4.  **`run_susie.R`**: The core fine-mapping script. It loads the GWAS
    summary statistics and pairs them with the PLINK2 correlation
    matrix, fits the `susie_rss` model, and exports 90% Credible Sets
    (CS) and Posterior Inclusion Probabilities (PIPs).

## Environment

    ## R version 4.1.2 (2021-11-01)
    ## Platform: x86_64-conda-linux-gnu (64-bit)
    ## Running under: Rocky Linux 8.10 (Green Obsidian)
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /programs/biogrids/x86_64-linux/r/4.1/lib/libopenblasp-r0.3.17.so
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] susieR_0.12.35    data.table_1.15.2 forcats_0.5.2     stringr_1.5.1    
    ##  [5] dplyr_1.1.4       purrr_1.2.2       readr_2.1.5       tidyr_1.3.1      
    ##  [9] tibble_3.2.1      ggplot2_3.5.1     tidyverse_1.3.2  
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.2.1    xfun_0.42           haven_2.5.1        
    ##  [4] gargle_1.5.2        lattice_0.20-45     colorspace_2.1-1   
    ##  [7] vctrs_0.6.5         generics_0.1.3      htmltools_0.5.7    
    ## [10] yaml_2.3.12         rlang_1.1.5         mixsqp_0.3-54      
    ## [13] pillar_1.10.2       glue_1.8.0          withr_3.0.2        
    ## [16] DBI_1.2.3.9005      dbplyr_2.4.0        modelr_0.1.9       
    ## [19] readxl_1.4.1        plyr_1.8.9          matrixStats_1.2.0  
    ## [22] lifecycle_1.0.4     munsell_0.5.1       gtable_0.3.6       
    ## [25] cellranger_1.1.0    rvest_1.0.3         evaluate_0.23      
    ## [28] knitr_1.45          tzdb_0.4.0          fastmap_1.1.1      
    ## [31] irlba_2.3.5.1       Rcpp_1.0.14         broom_1.0.0        
    ## [34] scales_1.3.0        backports_1.4.1     googlesheets4_1.0.1
    ## [37] jsonlite_1.8.8      fs_1.6.3            hms_1.1.3          
    ## [40] digest_0.6.35       stringi_1.8.3       grid_4.1.2         
    ## [43] cli_3.6.4           tools_4.1.2         magrittr_2.0.3     
    ## [46] crayon_1.5.2        pkgconfig_2.0.3     Matrix_1.6-3       
    ## [49] xml2_1.3.6          reprex_2.0.2        googledrive_2.0.0  
    ## [52] lubridate_1.8.0     reshape_0.8.9       rmarkdown_2.26     
    ## [55] httr_1.4.7          rstudioapi_0.15.0   R6_2.6.1           
    ## [58] compiler_4.1.2
