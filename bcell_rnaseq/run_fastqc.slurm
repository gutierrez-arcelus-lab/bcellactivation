#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=16gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --array=1-74
#SBATCH --job-name=FASTQC
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/FASTQC.%A-%a

source /programs/biogrids.shrc
export FASTQC_X=0.11.9

# INPUT
META=/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/B_cells_rnaseq/metadata.tsv
FQDIR=/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/B_cells_rnaseq/34.198.31.178/221025_MG10430_fastq
FQ1=$( awk 'FNR>1' $META | awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $4 }' )
FQ2=$( awk 'FNR>1' $META | awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $5 }' )

fastqc -t $SLURM_CPUS_PER_TASK ${FQDIR}/${FQ1}
fastqc -t $SLURM_CPUS_PER_TASK ${FQDIR}/${FQ2}

