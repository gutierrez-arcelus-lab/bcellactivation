#!/bin/bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=8gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --array=1-22
#SBATCH --job-name=vcf_merge
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/vcf_merge.%A.%a

# Load tools
source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

# Define chromosome
CHR=$SLURM_ARRAY_TASK_ID

# Define IO
OUTDIR=${SLURM_SUBMIT_DIR}/results/VCF

## VCF
VCFMGB=( $(echo ${OUTDIR}/chr${CHR}.MGB.04{01..08}.vcf.gz) )
VCFMGB+=( $(echo ${OUTDIR}/chr${CHR}.MGB.0410.vcf.gz) )
VCFKGP=$(echo ${OUTDIR}/chr${CHR}.1000G.vcf.gz)
ALLVCF=( "${VCFMGB[@]}" $VCFKGP)
NVCF=$(echo "${#ALLVCF[@]}" )
OUTINTER=${OUTDIR}/chr${CHR}_intersect

# Intersect MGB and 1000G VCFs for the same set of variants
bcftools isec -p $OUTINTER -n=$NVCF -c snps -O z "${ALLVCF[@]}"

# Merge files
VCFMERGED=${OUTDIR}/chr${CHR}.isec.merged.vcf.gz
VCFMAX=$(( $NVCF - 1 ))
INTERVCF=( $(eval echo ${OUTINTER}/000{0..$VCFMAX}.vcf.gz) )

bcftools merge --threads $SLURM_CPUS_PER_TASK "${INTERVCF[@]}" |\
    bcftools annotate -x 'INFO' -O z -o $VCFMERGED -

# Prune SNPs for r2 < 0.1
PRUNESNPS=${OUTDIR}/chr${CHR}.prunedsnps
OUTPREFIX=${OUTDIR}/chr${CHR}.isec.merged.pruned

plink --vcf $VCFMERGED \
    --indep-pairwise 50 10 0.1 \
    --out $PRUNESNPS

plink --vcf $VCFMERGED \
    --extract $PRUNESNPS.prune.in \
    --keep-allele-order \
    --recode vcf \
    --out $OUTPREFIX

