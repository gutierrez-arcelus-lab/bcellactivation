#!/bin/bash

#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --job-name=PCA
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/pca.%j

source /programs/biogrids.shrc

PREFIX=${SLURM_SUBMIT_DIR}/results/VCF/allchr.isec.merged.pruned
VCF=${PREFIX}.vcf.gz
OUT=${PREFIX}.pca

plink --vcf $VCF --pca --threads $SLURM_CPUS_PER_TASK --out $OUT
