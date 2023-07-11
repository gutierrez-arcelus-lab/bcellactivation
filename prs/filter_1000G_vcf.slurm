#!/bin/bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=8gb
#SBATCH -t 24:00:00
#SBATCH -p bch-compute
#SBATCH -a 1-22
#SBATCH -J bcftoolskg
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/bcftools_1kg-%A-%a.out
#SBATCH -e /temp_work/ch229163/log/bcftools_1kg-%A-%a.err

# Load tools
source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

# Define IO
CHR=chr${SLURM_ARRAY_TASK_ID}
VCFIN=/reference_databases/1000G_VCF/GRCh38/Genotype_VCFs/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_${CHR}.recalibrated_variants.vcf.gz
SAMPLES=${SLURM_SUBMIT_DIR}/mgb_data/kgp_samples.txt
OUTDIR=${TEMP_WORK}/vcf/prs
VCFOUT=${OUTDIR}/${CHR}.1000G.vcf.gz

bcftools view --threads $SLURM_CPUS_PER_TASK --samples-file $SAMPLES --force-samples $VCFIN |\
bcftools view --threads $SLURM_CPUS_PER_TASK --genotype ^miss --types snps -m2 -M2 --min-af 0.05:minor - |\
    bcftools annotate --threads $SLURM_CPUS_PER_TASK -x QUAL,FILTER,INFO,FORMAT --set-id +'%CHROM:%POS:%REF:%ALT' -O z -o $VCFOUT -
 
tabix -p vcf $VCFOUT
