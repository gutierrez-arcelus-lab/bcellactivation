#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=24gb
#SBATCH -t 24:00:00
#SBATCH -p bch-compute,bch-compute-pe
#SBATCH -J Salmon-index
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/Salmon-index.%j.out
#SBATCH -e /temp_work/ch229163/log/Salmon-index.%j.err

source /programs/biogrids.shrc
export SALMON_X=1.10.1

GENOME=/lab-share/IM-Gutierrez-e2/Public/References/Genomes/hsapiens/GRCh38/GRCh38.primary_assembly.genome.fa.gz
GENCODE=/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens/gencode.v41.primary_assembly.transcripts.fa.gz
GENTROME=${SLURM_SUBMIT_DIR}/data/salmon/gentrome.fa.gz 
DECOYS=${SLURM_SUBMIT_DIR}/data/salmon/decoys.txt
OUT=${SLURM_SUBMIT_DIR}/index/salmon

# Make decoys
grep "^>" <( gunzip -c $GENOME ) | cut -d " " -f 1 > $DECOYS
sed -i.bak -e 's/>//g' $DECOYS

# Make gentrome
cat $GENCODE $GENOME > $GENTROME

# Run indexing
salmon index \
    -t $GENTROME \
    -i $OUT \
    --decoys $DECOYS \
    -k 31 \
    -p $SLURM_CPUS_PER_TASK \
    --tmpdir ${TEMP_WORK}/temp/salmon
