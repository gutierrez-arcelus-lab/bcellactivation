#!/usr/bin/bash 

#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=16gb
#SBATCH -t 24:00:00
#SBATCH -p bch-compute,bch-compute-pe
#SBATCH -a 1-6
#SBATCH -J ldsc
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/ldsc.%A-%a.out
#SBATCH -e /temp_work/ch229163/log/ldsc.%A-%a.err

source /programs/biogrids.shrc
source $HOME/miniconda3/etc/profile.d/conda.sh

conda activate ldsc

cd $SLURM_SUBMIT_DIR

GWAS=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID' ./gwas_files.txt )
GWASID=$( basename $GWAS | sed 's/\.sumstats\.gz//' )

WEIGHTS=data/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.
FREQ=data/1000G_Phase3_frq/1000G.EUR.QC.

BASELINE=data/1000G_Phase3_baseline_v1.2_ldscores/baseline.
TURQ=data/gene_sets/ldscores/turquoise.
BLUE=data/gene_sets/ldscores/blue.
BROWN=data/gene_sets/ldscores/brown.
YELLOW=data/gene_sets/ldscores/yellow.
GREEN=data/gene_sets/ldscores/green.
RED=data/gene_sets/ldscores/red.
BLACK=data/gene_sets/ldscores/black.
CTS=data/gene_sets/module.ldcts

# Run S-LDSC
python2.7 /lab-share/IM-Gutierrez-e2/Public/tools/ldsc/ldsc.py \
    --h2 $GWAS \
    --ref-ld-chr $BASELINE,$TURQ,$BLUE,$BROWN,$YELLOW,$GREEN,$RED,$BLACK \
    --w-ld-chr $WEIGHTS \
    --frqfile-chr $FREQ \
    --overlap-annot --print-coefficients \
    --out results/modules_500_str_${GWASID}

# Run LDSC-SEG
python2.7 /lab-share/IM-Gutierrez-e2/Public/tools/ldsc/ldsc.py \
    --h2-cts $GWAS \
    --ref-ld-chr $BASELINE \
    --ref-ld-chr-cts $CTS \
    --w-ld-chr $WEIGHTS \
    --out results/modules_500_seg_${GWASID}

conda deactivate
