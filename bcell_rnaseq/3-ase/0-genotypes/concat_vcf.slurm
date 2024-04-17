#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=16gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute,bch-compute-pe
#SBATCH --array=1-16
#SBATCH --job-name=concat
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/concat-vcf-%A-%a.out
#SBATCH -e /temp_work/ch229163/log/concat-vcf-%A-%a.err

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

SPEC=${SLURM_SUBMIT_DIR}/data/vcf_samples.tsv
DONORID=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'NR==ARRID { print $1 }' $SPEC )
BATCH=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'NR==ARRID { print $2 }' $SPEC )

VCFS=( $(eval echo ${TEMP_WORK}/vcf/ase2/chr{{1..22},X}_${BATCH}_${DONORID}.vcf.gz) )
OUT=${SLURM_SUBMIT_DIR}/data/${DONORID}.vcf.gz

bcftools concat -O z -o $OUT "${VCFS[@]}"
tabix -p vcf $OUT
