#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=16gb
#SBATCH -t 24:00:00
#SBATCH -p bch-compute
#SBATCH -a 1-22
#SBATCH --job-name=vcf-merge
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/vcf-merge-%A-%a.out
#SBATCH -e /temp_work/ch229163/log/vcf-merge-%A-%a.err

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

# Params
CHR=chr${SLURM_ARRAY_TASK_ID}
CPUS=$SLURM_CPUS_PER_TASK
DIR=${TEMP_WORK}/vcf/langefeld 

# For each chr and each batch, save only the SNPs present in all of them
VCFS=( $(eval echo ${DIR}/${CHR}.MGB.04{03..07}.vcf.gz) )
N=$( echo "${#VCFS[@]}" )
OUTISEC=${DIR}/${CHR}_isec

bcftools isec --threads $CPUS -n=$N -c none -O z -p $OUTISEC "${VCFS[@]}"

# For each chr, merge all the batches for the same set of SNPs
MAX=$(( $N - 1 ))
VCFISEC=( $(eval echo ${OUTISEC}/000{0..$MAX}.vcf.gz) )
VCFMERGED=${DIR}/${CHR}.MGB.merged.vcf.gz

bcftools merge --threads $CPUS -O z -o $VCFMERGED "${VCFISEC[@]}"  

rm -r $OUTISEC
