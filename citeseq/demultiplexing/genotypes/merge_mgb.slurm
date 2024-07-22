#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=16gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute,bch-compute-pe
#SBATCH --array=1-23
#SBATCH --job-name=bcfmerge
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/bcftools_merge_%A_%a.out
#SBATCH -e /temp_work/ch229163/log/bcftools_merge_%A_%a.err

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

CHROM=$SLURM_ARRAY_TASK_ID

if [ "$CHROM" -eq 23 ]; then CHROM=X; fi

VCFS=( $(echo ${TEMP_WORK}/vcf/tmp/chr${CHROM}_04{03,{05..07}}.vcf.gz) )

N=$( echo "${#VCFS[@]}" )
OUTISEC=${TEMP_WORK}/vcf/tmp/chr${CHROM}_intersect

bcftools isec -n=$N -c none -O z -p $OUTISEC "${VCFS[@]}"

MAX=$(( $N - 1 ))
ISECVCFS=( $(eval echo ${OUTISEC}/000{0..$MAX}.vcf.gz) )
VCFMERGED=${TEMP_WORK}/vcf/tmp/chr${CHROM}_merged.vcf.gz

bcftools merge --threads $SLURM_CPUS_PER_TASK "${ISECVCFS[@]}" |\
    bcftools norm --threads $SLURM_CPUS_PER_TASK -m +both - |\
    bcftools view --threads $SLURM_CPUS_PER_TASK --types snps --genotype ^miss -m2 -M2 --min-ac 1:minor -e 'INFO/AC==INFO/AN' -O z -o $VCFMERGED - 
