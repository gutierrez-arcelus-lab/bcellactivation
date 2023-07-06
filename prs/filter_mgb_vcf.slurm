#!/bin/bash

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

# Define directory and chromosome given array spec file
SPEC=${SLURM_SUBMIT_DIR}/array.spec
BATCH=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $1 }' $SPEC )
CHR=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $2 }' $SPEC )

# Define IO
LABSHARE=/lab-share/IM-Gutierrez-e2/Public
VCFIN=${LABSHARE}/mgb_biobank/$BATCH/vcf/chr${CHR}.dose.vcf.gz 
IDS=${SLURM_SUBMIT_DIR}/mgb_data/ids_per_batch/${BATCH}.txt
OUTDIR=${TEMP_WORK}/vcf/prs
VCFOUT=${OUTDIR}/chr${CHR}.MGB.${BATCH}.vcf.gz

# CMD
bcftools view --threads $SLURM_CPUS_PER_TASK --samples-file $IDS --force-samples $VCFIN |\
    bcftools annotate --threads $SLURM_CPUS_PER_TASK -x INFO,FORMAT - |\
    bcftools view --threads $SLURM_CPUS_PER_TASK --min-ac 1:minor - |\
    bgzip -c > $VCFOUT    

tabix -p vcf $VCFOUT
