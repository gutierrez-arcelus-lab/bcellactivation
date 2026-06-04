# README

This repository contains the custom computational pipelines and R scripts used
for the data analysis and figure generation in the manuscript: *"A multi-omics
resource of human B-cell activation reveals dynamic response of
disease-associated loci"* (Aguiar, Franco, et al., 2026).

**Please navigate to each specific directory for detailed, step-by-step README instructions on how to run its respective pipeline.**

* `01_rnaseq_lowinput/` - Pipeline for quantification and analysis of the low-input RNA-seq timecourse data.
* `02_rnaseq_stdinput/` - Pipeline for the standard-input RNA-seq data, including differential expression and alternative splicing analyses.
* `03_atacseq/` - Upstream processing, peak calling, and analysis of chromatin accessibility (ATAC-seq) data.
* `04_citeseq/` - Single-cell multi-omics (CITE-seq) processing and analysis.
* `gwas_finemapping/` - Statistical fine-mapping of SLE GWAS loci using SuSiE.
* `paper_figures/` - R scripts to assemble, format, and export the final publication figures (e.g., Figures 1-6).
* `supp_data/` - Supplementary datasets generation. 

## 2. System Requirements

### Operating Systems

* **Upstream processing scripts (Bash/nf-core):** Executed on a Linux-based
  High-Performance Computing (HPC) environment.
* **Downstream analysis scripts (R):** Platform-independent. Tested on Linux.

### Software Dependencies

* **Software:**
    * Trim Galore v0.6.6
    * FastQC v0.11.9
    * PLINK v1.9
    * demuxlet "first release"
    * nf-core/atacseq v2.0
    * STAR v2.7.9a
    * Salmon v1.5.1
    * LDSC v1.0.1
    * HOMER v4.11
    * regtools v0.5.1
    * leafcutter v.0.2.9

* **R code and key packages:**
    * R v4.1.2
    * Key packages include:
	* `clusterProfiler` v4.4.2
	* `demuxmix` v1.6.0
	* `DESeq2` v1.34
	* `edgeR` v3.36
	* `WGCNA` v1.72
	* `Seurat` v4.1.1
	* `Harmony` v1.2.0
	* `SuSiE` v0.12.35
	* `fgsea` v1.20.0
	* `leafviz` v0.1.0

## 3. Installation Guide

* **Instructions:** This repository provides custom data analysis pipelines
  rather than a standalone software package. No formal installation is
  required. Users must independently install R and the third-party tools listed
  above according to their official documentation.

## 4. Demo

* **Instructions & Expected Output:** Because this repository contains custom
  analysis scripts tailored to a specific project rather than generalized
  software, a simulated demo dataset is not provided. To reproduce the
  analysis, users should download the raw and processed datasets associated
  with this manuscript from dbGaP and Zenodo and execute the provided scripts.

* **Expected Run Time:** Upstream steps such as raw read processing and
  alignment run times depend heavily on HPC resources and data size (typically
  several hours to days). Downstream R scripts typically run in 5 to 30 minutes
  on a standard desktop computer.

## 5. Instructions for Use

1.  **Download Data:** Obtain the corresponding datasets from dbGaP
(phs004221.v1.p1) and Zenodo (https://zenodo.org/uploads/15595470)
2.  **Execute:** Run the scripts sequentially or as needed to reproduce
specific steps of the analysis.
