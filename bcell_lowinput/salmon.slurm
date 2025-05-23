#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 6
#SBATCH --mem=24gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --array=1-174
#SBATCH --job-name=Salmon
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/Salmon-low.%A-%a

source /programs/biogrids.shrc
export SALMON_X=1.10.1

# INPUT
SAMPLELIST=/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/B_cells_lowinput/metadata/description_hsapiens_qc.tsv
PLATE=$( awk 'FNR>1' $SAMPLELIST | awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $1 }')
WELL=$( awk 'FNR>1' $SAMPLELIST | awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $2 }')
FQ1=$( awk 'FNR>1' $SAMPLELIST | awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $7 }')
FQ2=$( awk 'FNR>1' $SAMPLELIST | awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $8 }')
FQLIST1=( $(echo $FQ1 | tr "," " ") )
FQLIST1=( "${FQLIST1[@]/#/${TEMP_WORK}/fastq/lowinput/}" )
FQLIST2=( $(echo $FQ2 | tr "," " ") )
FQLIST2=( "${FQLIST2[@]/#/${TEMP_WORK}/fastq/lowinput/}" )

# Index
INDEX=${HOME}/lupus/indices/SALMON

# Output
OUT=${HOME}/lupus/bcell_lowinput/results/salmon/${PLATE}_${WELL}
mkdir -p $OUT

# Run Salmon quantification
salmon quant -i $INDEX -l A \
    -1 ${FQLIST1[@]} -2 ${FQLIST2[@]} \
    -p $SLURM_CPUS_PER_TASK \
    -o $OUT \
    --seqBias --posBias --gcBias
