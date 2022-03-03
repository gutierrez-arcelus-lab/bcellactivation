#!/bin/bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=4gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --array=0-22
#SBATCH --job-name=bcftools-stats
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/bcftools_stats.%A.%a

# Load tools
source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

CHRS=( $(seq 1 22; echo X) )
CHR=${CHRS["$SLURM_ARRAY_TASK_ID"]}

VCF=/lab-share/IM-Gutierrez-e2/Public/mgb_biobank/0410/vcf/chr${CHR}.dose.vcf.gz
OUT=${TEMP_WORK}/VCF/chr${CHR}.MGB.0410.stats

bcftools stats --threads $SLURM_CPUS_PER_TASK $VCF > $OUT
