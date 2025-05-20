#!/bin/bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=8gb
#SBATCH -t 24:00:00
#SBATCH -p bch-compute
#SBATCH -J PCA-1KG
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/pca_1kg-%A-%a.out
#SBATCH -e /temp_work/ch229163/log/pca_1kg-%A-%a.err

# Load tools
source /programs/biogrids.shrc
export BCFTOOLS_X=1.12
export PLINK2_X=2.00a3

# Define IO
VCFIN=${TEMP_WORK}/vcf/prs/allchr.merged.pruned.vcf.gz
SAMPLES=${SLURM_SUBMIT_DIR}/data/kgp_ids.txt
OUTDIR=${TEMP_WORK}/vcf/prs
VCFOUT=${OUTDIR}/allchr.1000G.vcf.gz
PCAOUT=${OUTDIR}/allchr.1000G.plink.pca

bcftools view --threads $SLURM_CPUS_PER_TASK --samples-file $SAMPLES --force-samples $VCFIN |\
    bcftools view --threads $SLURM_CPUS_PER_TASK --min-af 0.05:minor -O z -o $VCFOUT -
 
tabix -p vcf $VCFOUT

plink2 --vcf $VCFOUT --pca approx --threads $SLURM_CPUS_PER_TASK --out $PCAOUT
