#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=16gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --array=1-115
#SBATCH --job-name=vcf-sub
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/vcf-subset-%A-%a.out
#SBATCH -e /temp_work/ch229163/log/vcf-subset-%A-%a.err

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

SPEC=/lab-share/IM-Gutierrez-e2/Public/vitor/lupus/bcell_rnaseq/4-fix_ase/0-genotypes/array_spec.txt
BATCH=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'NR==ARRID { print $1 }' $SPEC )
CHR=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'NR==ARRID { print $2 }' $SPEC )

VCFIN=/lab-share/IM-Gutierrez-e2/Public/External_datasets/mgb_biobank/${BATCH}/vcf/${CHR}.dose.vcf.gz 
IDS=/lab-share/IM-Gutierrez-e2/Public/vitor/lupus/mgb_biobank/results/eur_females_${BATCH}.txt
VCFVAR=${SLURM_SUBMIT_DIR}/data/variants.vcf.gz
PREFIX=${TEMP_WORK}/vcf/quasar/${CHR}_${BATCH} 
VCFOUT=${PREFIX}.vcf.gz
     
bcftools isec --threads $SLURM_CPUS_PER_TASK -c none -n=2 -w2 $VCFVAR $VCFIN |\
    bcftools view --threads $SLURM_CPUS_PER_TASK --samples-file $IDS --force-samples - |\
    bcftools annotate --threads $SLURM_CPUS_PER_TASK -x INFO,FORMAT -O z -o $VCFOUT -

tabix -p vcf $VCFOUT
