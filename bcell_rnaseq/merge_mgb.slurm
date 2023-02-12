#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=16gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --array=22
#SBATCH --job-name=mergeVCF
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/mergeVCF-%A-%a

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

CHRS=( $(seq 1 22; echo X) )
CHR=chr${CHRS[$SLURM_ARRAY_TASK_ID]}
CPUS=$SLURM_CPUS_PER_TASK

VCFS=( $(eval echo ${TEMP_WORK}/${CHR}.MGB.04{03..07}.vcf.gz) )
N=$( echo "${#VCFS[@]}" )
OUTISEC=${TEMP_WORK}/${CHR}_intersect

bcftools isec --threads $CPUS -n=$N -c none -O z -p $OUTISEC "${VCFS[@]}"

MAX=$(( $N - 1 ))
ISECVCFS=( $(eval echo ${OUTISEC}/000{0..$MAX}.vcf.gz) )
VCFMERGED=${TEMP_WORK}/${CHR}.MGB.merged.vcf.gz

bcftools merge --threads $CPUS "${ISECVCFS[@]}" |\
    bcftools annotate --threads $CPUS -x INFO,FORMAT - |\
    bcftools view --threads $CPUS --min-ac 1:minor -O z -o $VCFMERGED - 
