#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=16gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --array=1-22
#SBATCH --job-name=GetMGB
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/mgb-%A-%a

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

CHR=$SLURM_ARRAY_TASK_ID
INFO=${SLURM_SUBMIT_DIR}/samples.tsv
BATCHES=$( awk '{ print $1 }' $INFO | sort | uniq )

for b in $BATCHES
do
    VCFIN=/lab-share/IM-Gutierrez-e2/Public/External_datasets/mgb_biobank/${b}/vcf/chr${CHR}.dose.vcf.gz 
    IDS=${SLURM_SUBMIT_DIR}/data/${b}.txt
    VCFOUT=${TEMP_WORK}/chr${CHR}.MGB.${b}.vcf.gz

    bcftools view --threads $SLURM_CPUS_PER_TASK --samples-file $IDS --force-samples $VCFIN |\
	bcftools norm --threads $SLURM_CPUS_PER_TASK -m +both - |\
	bcftools view --threads $SLURM_CPUS_PER_TASK --genotype ^miss --types snps -m2 -M2 -O z -o $VCFOUT -

    tabix -f -p vcf $VCFOUT

    VCFS+=( $VCFOUT )
done

N=$(echo "${#VCFS[@]}" )
OUTISEC=${TEMP_WORK}/chr${CHR}_intersect

bcftools isec -n=$N -c none -O z -p $OUTISEC "${VCFS[@]}"

MAX=$(( $N - 1 ))
ISECVCFS=( $(eval echo ${OUTISEC}/000{0..$MAX}.vcf.gz) )
VCFMERGED=${TEMP_WORK}/chr${CHR}.MGB.merged.vcf.gz

bcftools merge --threads $SLURM_CPUS_PER_TASK "${ISECVCFS[@]}" |\
    bcftools view --threads $SLURM_CPUS_PER_TASK --min-ac 1:minor -O z -o $VCFMERGED - 
