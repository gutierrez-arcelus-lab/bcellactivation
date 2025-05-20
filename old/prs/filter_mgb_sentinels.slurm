#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=8gb
#SBATCH -t 24:00:00
#SBATCH -p bch-compute
#SBATCH -a 1-207
#SBATCH -J bcftools
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/bcftools-%A-%a.out
#SBATCH -e /temp_work/ch229163/log/bcftools-%A-%a.err

# Load tools
source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

cd $SLURM_SUBMIT_DIR

# Define batch and chromosome given array spec file
SPEC=./array.spec
BATCH=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $1 }' $SPEC )
CHR=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $2 }' $SPEC )

# Define IO
VCFDIR=/lab-share/IM-Gutierrez-e2/Public/External_datasets/mgb_biobank/${BATCH}/vcf
VCFIN=${VCFDIR}/${CHR}.dose.vcf.gz 
IDS=./mgb_data/ids_per_batch/${BATCH}.txt
SNPS=./ld_partners/hg19/data/sentinels_grch38.tsv
OUTDIR=${TEMP_WORK}/vcf/prs
VCFOUT=${OUTDIR}/${CHR}.MGB.${BATCH}.sentinels.vcf.gz

# CMD
bcftools view --threads $SLURM_CPUS_PER_TASK -R $SNPS $VCFIN |\
    bcftools view --samples-file $IDS --force-samples - |\
    bcftools view --threads $SLURM_CPUS_PER_TASK -i 'INFO/R2>=0.8' --min-ac 1:minor - |\
    bcftools annotate --threads $SLURM_CPUS_PER_TASK -x FILTER,INFO,FORMAT -O z -o $VCFOUT -

tabix -p vcf $VCFOUT
