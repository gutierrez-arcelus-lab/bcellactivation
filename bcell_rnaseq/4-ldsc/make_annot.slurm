#!/usr/bin/bash 

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=8gb
#SBATCH -t 24:00:00
#SBATCH -p bch-compute
#SBATCH -a 1-110
#SBATCH -J ldsc
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/ldsc-annot.%A-%a.out
#SBATCH -e /temp_work/ch229163/log/ldsc-annot.%A-%a.err

source /programs/biogrids.shrc
source $HOME/miniconda3/etc/profile.d/conda.sh

conda activate ldsc

cd $SLURM_SUBMIT_DIR

# Input/Output
SPEC=data/array_spec.tsv
CHROM=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $1 }' $SPEC)
SET=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $2 }' $SPEC)

# LDSC data
BIM=/lab-share/IM-Gutierrez-e2/Public/tools/ldsc_data/v4/1000G_EUR_Phase3_plink/1000G.EUR.QC.${CHROM}.bim
HAPMAP=/lab-share/IM-Gutierrez-e2/Public/tools/ldsc_data/v4/hm3_no_MHC.list.txt

# Gene/Peak set data
COORD=data/gene_coord.txt
GENESET=data/${SET}.txt
OUT=data/ldscores/${SET}.${CHROM}
ANNOT=${OUT}.annot

# Make annotation
python.pybedtools /lab-share/IM-Gutierrez-e2/Public/tools/ldsc/make_annot.py \
    --gene-set-file $GENESET \
    --gene-coord-file $COORD \
    --bimfile $BIM \
    --windowsize 100000 \
    --annot-file $ANNOT

# Compute LD scores
python2.7 /lab-share/IM-Gutierrez-e2/Public/tools/ldsc/ldsc.py \
    --l2 \
    --bfile ${BIM%.*} \
    --annot $ANNOT \
    --print-snps $HAPMAP \
    --ld-wind-cm 1 \
    --thin-annot \
    --out $OUT

# Rename annotation
zcat ${OUT}.l2.ldscore.gz | sed "1,1s/L2/L2${SET}/" | gzip > ${OUT}.l2.tmp.gz
mv ${OUT}.l2.tmp.gz ${OUT}.l2.ldscore.gz 

conda deactivate
