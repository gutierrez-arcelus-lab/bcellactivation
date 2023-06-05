#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=16gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --array=0-22
#SBATCH --job-name=mergeVCF
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/mergeVCF-%A-%a

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

CHRS=( $(seq 1 22; echo X) )
CHR=chr${CHRS[$SLURM_ARRAY_TASK_ID]}
CPUS=$SLURM_CPUS_PER_TASK

# For each chr and each batch, save only the SNPs present in all of them
VCFS=( $(eval echo ${TEMP_WORK}/vcf/mbv/${CHR}.MGB.040{3,4,6,7}.vcf.gz) )
N=$( echo "${#VCFS[@]}" )
OUTISEC=${TEMP_WORK}/vcf/mbv/${CHR}_intersect

bcftools isec --threads $CPUS -n=$N -c none -O z -p $OUTISEC "${VCFS[@]}"

# For each chr, merge all the batches for the same set of SNPs
MAX=$(( $N - 1 ))
ISECVCFS=( $(eval echo ${OUTISEC}/000{0..$MAX}.vcf.gz) )
VCFMERGED=${TEMP_WORK}/vcf/mbv/${CHR}.MGB.merged.vcf.gz

bcftools merge --threads $CPUS "${ISECVCFS[@]}" |\
    bcftools annotate --threads $CPUS -x INFO,FORMAT - |\
    bcftools view --threads $CPUS --min-ac 1:minor -O z -o $VCFMERGED - 
