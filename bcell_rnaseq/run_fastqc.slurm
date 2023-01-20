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
export CUTADAPT_X=4.1
export FASTQC_X=0.11.9
export PYTHON_X=3.8.8

# INPUT
META=/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/B_cells_rnaseq/metadata.tsv
FQDIR=/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/B_cells_rnaseq/34.198.31.178/221025_MG10430_fastq
FQ1=$( awk 'FNR>1' $META | awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $4 }' )
FQ2=$( awk 'FNR>1' $META | awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $5 }' )

#fastqc -t $SLURM_CPUS_PER_TASK ${FQDIR}/${FQ1}
#fastqc -t $SLURM_CPUS_PER_TASK ${FQDIR}/${FQ2}

QCDIR=/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/B_cells_rnaseq/QC
mkdir -p $QCDIR

trim_galore --paired --length 50 --illumina -j 4 -o $QCDIR ${FQDIR}/${FQ1} ${FQDIR}/${FQ2}

FQ1qc=$( echo $FQ1 | sed 's|.fastq.gz|_val_1.fq.gz|' )
FQ2qc=$( echo $FQ2 | sed 's|.fastq.gz|_val_2.fq.gz|' )

fastqc -t $SLURM_CPUS_PER_TASK ${QCDIR}/${FQ1qc}
fastqc -t $SLURM_CPUS_PER_TASK ${QCDIR}/${FQ2qc}

