CITE-seq Pilots
================

Raw counts
----------

![](README_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

Distribution of HTO counts
--------------------------

### Pilot 1

![](README_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

![](README_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

QC
--

Here we use the miQC package to model the percentage of mitochondrial
reads and number of genes, in order to identify and remove compromised
cells.

![](README_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

    # Removing 4104 out of 13946 cells.

    # Removing 723 out of 10562 cells.

Demultiplex cells based on HTO
------------------------------

![](README_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

![](README_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

![](README_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

In the plot above, we see 2 differences between pilots 1 and 2:

-   In pilot 1, HTO counts in higher, including for those cells being
    classified as “Negatives”;

-   In pilot 2, HTO counts for doublets not that higher than those for
    singlets.

Regarding the last point, doublets don’t necessarily have higher counts,
but a mixture of HTOs.

However, a lot of droplets being classified as doublets in pilot 2 have
a profile resembling singlets in respect to the mixture of HTO, as we
can see in the plots below.

![](README_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->
![](README_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

Demultiplex by individuals’ genotypes or GMM-demux in Pilot 2
-------------------------------------------------------------

![](README_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

Scrublet
--------

![](README_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

Custom demultiplexing based on HTO
----------------------------------

Number of cells for custom vs Seurat demultiplex method
-------------------------------------------------------

![](README_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

Profile of Singlets and Doublets in custom demultiplexing
---------------------------------------------------------

![](README_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

![](README_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

![](README_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

Integrate Pilot 1 and Pilot 2
-----------------------------

![](README_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->
