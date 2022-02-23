#!/bin/bash

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=8gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --array=1-22
#SBATCH --job-name=filter_1000G
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/filter_1000G.%A.%a

# Load tools
source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

# Define IO
CHR=$SLURM_ARRAY_TASK_ID
LABSHARE=/lab-share/IM-Gutierrez-e2/Public
KGDIR=${LABSHARE}/vcf_1000G
VCFIN=${KGDIR}/ALL.chr${CHR}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz

OUTDIR=${TEMP_WORK}/VCF
VCFTMP=${OUTDIR}/chr${CHR}.1000G.tmp.vcf.gz
VCFOUT=${OUTDIR}/chr${CHR}.1000G.vcf.gz
VARSKG=${OUTDIR}/chr${CHR}.1000G.variants.txt

# first adjust chr notation
zcat $VCFIN |\
    awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' - |\
    bgzip -c > $VCFTMP

tabix -p vcf $VCFTMP

# Filter 1000G
bcftools view --genotype ^miss --types snps -m2 -M2 $VCFTMP |\
    bcftools view --min-af 0.01:minor - |\
    bcftools annotate -x INFO,FORMAT - |\
    bgzip -c > $VCFOUT

tabix -p vcf $VCFOUT

# Save Chromosome, Position, ID, Ref and Alt alleles to a file
zcat $VCFOUT |\
    grep -v '^##' |\
    awk '{print $1, $2, $3, $4, $5}' > $VARSKG

