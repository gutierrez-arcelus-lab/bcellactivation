#!/usr/bin/sh

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=16gb
#SBATCH -t 24:00:00
#SBATCH -p bch-compute
#SBATCH -a 1-450
#SBATCH -J LD1KGP
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/LD1KGP-%A-%a.out
#SBATCH -e /temp_work/ch229163/log/LD1KGP-%A-%a.err

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

SAMPLEIDS=${SLURM_SUBMIT_DIR}/data/kgp_samples.txt
REGIONSFILE=${SLURM_SUBMIT_DIR}/data/regions.tsv 
REGION=$( awk -v ARR="$SLURM_ARRAY_TASK_ID" 'FNR == ARR { print $1 }' $REGIONSFILE )
RSID=$( awk -v ARR="$SLURM_ARRAY_TASK_ID" 'FNR == ARR { print $2 }' $REGIONSFILE )
CHR=$( echo $REGION | cut -d':' -f1 )
VCF1K=/reference_databases/1000G_VCF/phase3/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
PREFIX=${TEMP_WORK}/vcf/prs/hg19/${RSID}
VCFOUT=${PREFIX}.vcf.gz

bcftools view --threads $SLURM_CPUS_PER_TASK -r $REGION $VCF1K |\
    bcftools view --threads $SLURM_CPUS_PER_TASK --samples-file $SAMPLEIDS --force-samples - |\
    bcftools view --threads $SLURM_CPUS_PER_TASK --genotype ^miss -e INFO/AC==INFO/AN --min-ac 1:minor - |\
    bcftools annotate --threads $SLURM_CPUS_PER_TASK -x INFO,FORMAT -O z -o $VCFOUT -

tabix -p vcf $VCFOUT

GWASPOS=${SLURM_SUBMIT_DIR}/data/gwas_positions.tsv
VCFFILT=${PREFIX}.filt.vcf.gz

bcftools view --threads $SLURM_CPUS_PER_TASK -R $GWASPOS -O z -o $VCFFILT $VCFOUT

plink --vcf $VCFFILT --r2 square spaces --keep-allele-order --out ${PREFIX}
rm ${PREFIX}.nosex
