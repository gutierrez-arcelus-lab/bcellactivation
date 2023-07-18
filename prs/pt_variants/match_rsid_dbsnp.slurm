#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=16gb
#SBATCH -t 24:00:00
#SBATCH -p bch-compute
#SBATCH -J bcftools
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/bcftools-dbSNP-%A-%a.out
#SBATCH -e /temp_work/ch229163/log/bcftools-dbSNP-%A-%a.err

# Load tools
source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

# Define IO
DB=/lab-share/IM-Gutierrez-e2/Public/References/dbSNP/GCF_000001405.39.gz 
OUT=${SLURM_SUBMIT_DIR}/data/pt_5e-4.vcf

# CMD
bcftools annotate --threads $SLURM_CPUS_PER_TASK -x INFO $DB |\
    bcftools view --threads $SLURM_CPUS_PER_TASK -i "%ID=@/lab-share/IM-Gutierrez-e2/Public/vitor/lupus/prs/pt_variants/data/pt_5e-4_rsids.txt" -O v -o $OUT -

