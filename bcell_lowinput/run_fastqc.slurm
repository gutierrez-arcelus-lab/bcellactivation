#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=8gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --array=1-732
#SBATCH --job-name=FASTQC
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/FASTQC.%A-%a

source /programs/biogrids.shrc
export FASTQC_X=0.11.9

# INPUT
SAMPLELIST=/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/B_cells_lowinput/metadata/description_hsapiens_longformat.tsv
FQ1=$( awk 'FNR>1' $SAMPLELIST | awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $8 }')
FQ2=$( awk 'FNR>1' $SAMPLELIST | awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $9 }')

fastqc -t $SLURM_CPUS_PER_TASK $FQ1
fastqc -t $SLURM_CPUS_PER_TASK $FQ2

