#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=16gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --array=1-368
#SBATCH --job-name=bcftools
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/subset-mgb-%A-%a.out
#SBATCH -e /temp_work/ch229163/log/subset-mgb-%A-%a.err

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

SPEC=${SLURM_SUBMIT_DIR}/array_spec.txt
DONORID=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'NR==ARRID { print $1 }' $SPEC )
SAMPLEID=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'NR==ARRID { print $2 }' $SPEC )
BATCH=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'NR==ARRID { print $3 }' $SPEC )
CHROM=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'NR==ARRID { print $4 }' $SPEC )

VCFIN=/lab-share/IM-Gutierrez-e2/Public/External_datasets/mgb_biobank/${BATCH}/vcf/${CHROM}.dose.vcf.gz 
VCFOUT=${TEMP_WORK}/vcf/ase2/${CHROM}_${DONORID}.vcf.gz

bcftools view --threads $SLURM_CPUS_PER_TASK --samples $SAMPLEID $VCFIN |\
    bcftools norm --threads $SLURM_CPUS_PER_TASK -m +both - |\
    bcftools view --threads $SLURM_CPUS_PER_TASK -i 'INFO/R2>=0.9' --types snps --genotype het -m2 -M2 - |\
    bcftools annotate --threads $SLURM_CPUS_PER_TASK -x FILTER,INFO,FORMAT -O z -o $VCFOUT -

tabix -p vcf $VCFOUT
