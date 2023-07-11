#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=16gb
#SBATCH -t 24:00:00
#SBATCH -p bch-compute
#SBATCH -a 1-22
#SBATCH -J vcf_merge
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/vcf_merge-%A-%a.out
#SBATCH -e /temp_work/ch229163/log/vcf_merge-%A-%a.err

# Load tools
source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

# Define chromosome
CHR=chr${SLURM_ARRAY_TASK_ID}

# Define IO
OUTDIR=${TEMP_WORK}/vcf/prs
VCFMGB=( $(echo ${OUTDIR}/${CHR}.MGB.04{{01..08},10}.vcf.gz) )
VCFKGP=$(echo ${OUTDIR}/${CHR}.1000G.vcf.gz)
ALLVCF=( "${VCFMGB[@]}" $VCFKGP)
NVCF=$(echo "${#ALLVCF[@]}" )
OUTISEC=${OUTDIR}/${CHR}_intersect

# Intersect MGB and 1000G VCFs for the same set of variants
bcftools isec -n=$NVCF -c none -O z -p $OUTISEC "${ALLVCF[@]}"

# Merge files
VCFMERGED=${OUTDIR}/${CHR}.merged.vcf.gz
VCFMAX=$(( $NVCF - 1 ))
INTERVCF=( $(eval echo ${OUTISEC}/000{0..$VCFMAX}.vcf.gz) )

bcftools merge --threads $SLURM_CPUS_PER_TASK "${INTERVCF[@]}" |\
    bcftools annotate -x 'INFO' -O z -o $VCFMERGED -

# Prune SNPs for r2 < 0.1
PRUNESNPS=${OUTDIR}/${CHR}.prunedsnps
OUTPREFIX=${OUTDIR}/${CHR}.merged.pruned

plink --vcf $VCFMERGED \
    --indep-pairwise 100kb 1 0.1 \
    --threads $SLURM_CPUS_PER_TASK \
    --out $PRUNESNPS

plink --vcf $VCFMERGED \
    --extract $PRUNESNPS.prune.in \
    --real-ref-alleles \
    --recode vcf \
    --threads $SLURM_CPUS_PER_TASK \
    --out $OUTPREFIX

