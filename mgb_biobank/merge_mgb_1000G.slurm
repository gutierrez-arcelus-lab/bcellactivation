#!/bin/bash

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=8gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --array=1-22
#SBATCH --job-name=vcf_merge
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /home/ch229163/log/vcf_merge.%A.%a

# Load tools
source /programs/biogrids.shrc

export BCFTOOLS_X=1.12
export R_X=4.1

# Define chromosome
CHR=$SLURM_ARRAY_TASK_ID

# Define IO
OUTDIR=$TEMP_WORK/VCF
VCFMGB=$OUTDIR/chr$CHR.MGB.04{01..04}.vcf.gz
VCFKGP=$OUTDIR/chr$CHR.1000G.vcf.gz
VARSMGB=$OUTDIR/chr$CHR.MGB.04{01..04}.variants.txt
VARS1KG=$OUTDIR/chr$CHR.1000G.variants.txt
OUTPOS=$OUTDIR/chr$CHR.commonvars.pos

# Extract common variants among datasets
Rscript $SLURM_SUBMIT_DIR/match_snps.R $( eval echo $VARSMGB ) $VARS1KG $OUTPOS

# Filter MGB and 1000G VCFs for the same set of variants
for f in $( eval echo $VCFMGB; echo $VCFKGP)
do
    VCFOUT=$( echo $f | sed 's/vcf\.gz/filtered.vcf.gz/' )
    VCFOUTLIST+=($VCFOUT)
    bcftools view -R $OUTPOS $f -O z -o $VCFOUT
    tabix -p vcf $VCFOUT
done

# Merge files
VCFMERGED=$OUTDIR/chr${CHR}.merged.vcf.gz

bcftools merge "${VCFOUTLIST[@]}" |\
    bcftools annotate -x 'INFO' - -O z -o $VCFMERGED

# Prune SNPs for r2 < 0.1
PRUNESNPS=$OUTDIR/chr${CHR}.prunedsnps
OUTPREFIX=$OUTDIR/chr${CHR}.merged.pruned

plink --vcf $VCFMERGED \
    --indep-pairwise 50 10 0.1 \
    --out $PRUNESNPS

plink --vcf $VCFMERGED \
    --extract $PRUNESNPS.prune.in \
    --keep-allele-order \
    --recode vcf \
    --out $OUTPREFIX

