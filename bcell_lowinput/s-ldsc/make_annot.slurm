#!/usr/bin/bash 

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=8gb
#SBATCH -t 24:00:00
#SBATCH -p bch-compute
#SBATCH -a 1-153
#SBATCH -J ldsc-annot
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/ldsc-annot.%A-%a.out
#SBATCH -e /temp_work/ch229163/log/ldsc-annot.%A-%a.err

source /programs/biogrids.shrc

cd $SLURM_SUBMIT_DIR

SPEC=data/gene_sets/array_spec.tsv
CHROM=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $1 }' $SPEC)
MODULE=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $2 }' $SPEC)

BIM=data/1000G_Phase3_plinkfiles/1000G_EUR_Phase3_plink/1000G.EUR.QC.${CHROM}.bim
COORD=data/gene_sets/ENSG_coord.txt
GENESET=data/gene_sets/genes/${MODULE}.txt
OUT=data/gene_sets/annot/${MODULE}.${CHROM}.annot.gz

python.pybedtools /lab-share/IM-Gutierrez-e2/Public/vitor/software/ldsc/make_annot.py \
    --gene-set-file $GENESET \
    --gene-coord-file $COORD \
    --bimfile $BIM \
    --windowsize 100000 \
    --annot-file $OUT
