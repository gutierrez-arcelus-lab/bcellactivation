#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=32gb
#SBATCH --time=24:00:00 
#SBATCH -p bch-compute 
#SBATCH --array=1-74
#SBATCH --job-name=Salmon
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/Salmon-quant-%A-%a

source /programs/biogrids/shrc.shrc
export SALMON_X=1.5.1

# INPUT
FQDIR=/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/B_cells_rnaseq/QC 
META=/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/B_cells_rnaseq/metadata.tsv
SAMPLEID=$( awk 'FNR>1' $META | awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'NR==ARRID { print $2 }' )
STIM=$( awk 'FNR>1' $META | awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'NR==ARRID { print $3 }' )
FQ1=$( awk 'FNR>1' $META | awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'NR==ARRID { print $4 }' )
FQ2=$( awk 'FNR>1' $META | awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'NR==ARRID { print $5 }' )
FQ1qc=${FQDIR}/$( echo $FQ1 | sed 's|.fastq.gz|_val_1.fq.gz|' )
FQ2qc=${FQDIR}/$( echo $FQ2 | sed 's|.fastq.gz|_val_2.fq.gz|' )

# INDEX
INDEX=/lab-share/IM-Gutierrez-e2/Public/vitor/lupus/indices/SALMON

# OUTPUT
OUT=${SLURM_SUBMIT_DIR}/results/${SAMPLEID}_${STIM}
mkdir -p $OUT

# EXECUTE
salmon quant -i $INDEX -l A \
    -1 $FQ1qc -2 $FQ2qc \
    -p $SLURM_CPUS_PER_TASK \
    -o $OUT \
    --seqBias --posBias --gcBias

