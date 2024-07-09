#!/usr/bin/bash 

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=8gb
#SBATCH -t 24:00:00
#SBATCH -p bch-compute,bch-compute-pe
#SBATCH -a 1-176
#SBATCH -J ldsc
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/ldsc-annot.%A-%a.out
#SBATCH -e /temp_work/ch229163/log/ldsc-annot.%A-%a.err

source /programs/biogrids.shrc
source $HOME/miniconda3/etc/profile.d/conda.sh

conda activate ldsc

cd $SLURM_SUBMIT_DIR

SPEC=data/gene_sets/array_spec.tsv
CHROM=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $1 }' $SPEC)
MODULE=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $2 }' $SPEC)

BIM=data/1000G_Phase3_plinkfiles/1000G_EUR_Phase3_plink/1000G.EUR.QC.${CHROM}.bim
COORD=data/gene_sets/ENSG_coord.txt
GENESET=data/gene_sets/genes/${MODULE}.txt
HAPMAP=data/print_snps.txt
OUT=data/gene_sets/ldscores/${MODULE}.${CHROM}
ANNOT=${OUT}.annot

# Make annotation
python.pybedtools /lab-share/IM-Gutierrez-e2/Public/vitor/software/ldsc/make_annot.py \
    --gene-set-file $GENESET \
    --gene-coord-file $COORD \
    --bimfile $BIM \
    --windowsize 100000 \
    --annot-file $ANNOT

# Compute LD scores
python2.7 /lab-share/IM-Gutierrez-e2/Public/vitor/software/ldsc/ldsc.py \
    --l2 \
    --bfile ${BIM%.*} \
    --annot $ANNOT \
    --print-snps $HAPMAP \
    --ld-wind-cm 1 \
    --thin-annot \
    --out $OUT

# Rename annotation
zcat ${OUT}.l2.ldscore.gz | sed "1,1s/L2/L2${MODULE}/" | gzip > ${OUT}.l2.tmp.gz
mv ${OUT}.l2.tmp.gz ${OUT}.l2.ldscore.gz 

conda deactivate
