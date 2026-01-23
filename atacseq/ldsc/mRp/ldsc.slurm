#!/usr/bin/bash 

#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=8gb
#SBATCH -t 01:00:00
#SBATCH -p bch-compute
#SBATCH -a 1-37
#SBATCH -J ldsc
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o logs/ldsc.%A-%a.out
#SBATCH -e logs/ldsc.%A-%a.err

source /programs/biogrids.shrc

cd $SLURM_SUBMIT_DIR

METADATA=./data/traits.txt
STATSDIR=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" -F '\t' 'FNR==ARRID { print $1 }' $METADATA ) 
TRAIT=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" -F '\t' 'FNR==ARRID { print $3 }' $METADATA )

LDSCDIR=/lab-share/IM-Gutierrez-e2/Public/tools/ldsc_data/v4 
GWAS=${LDSCDIR}/${STATSDIR}/${TRAIT}.sumstats.gz 
WEIGHTS=${LDSCDIR}/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.
FREQ=${LDSCDIR}/1000G_Phase3_frq/1000G.EUR.QC. 
BASELINE=${LDSCDIR}/1000G_Phase3_baseline_v1.2_ldscores/baseline.

IL4=./data/ldscores/IL4.
TLR=./data/ldscores/TLR7.
BCR=./data/ldscores/BCR.
DN2=./data/ldscores/DN2.
CTS=./data/peaksets.ldcts

# Run S-LDSC
${HOME}/miniconda3/envs/ldsc/bin/python /lab-share/IM-Gutierrez-e2/Public/tools/ldsc/ldsc.py \
    --h2 $GWAS \
    --ref-ld-chr $BASELINE,$IL4,$TLR,$BCR,$DN2 \
    --w-ld-chr $WEIGHTS \
    --frqfile-chr $FREQ \
    --overlap-annot --print-coefficients \
    --out ./results/sldsc_${TRAIT}

# Run LDSC-SEG
${HOME}/miniconda3/envs/ldsc/bin/python /lab-share/IM-Gutierrez-e2/Public/tools/ldsc/ldsc.py \
    --h2-cts $GWAS \
    --ref-ld-chr $BASELINE \
    --ref-ld-chr-cts $CTS \
    --w-ld-chr $WEIGHTS \
    --out ./results/ldscseg_${TRAIT}
