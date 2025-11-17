#!/usr/bin/env bash

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12
export GATK_X=4.1.4.1

LABSHR=/lab-share/IM-Gutierrez-e2/Public
ID=MG8989
OUT=${LABSHR}/vitor/ase/bcell_bulk/results/gatk/${ID}.knownSNVs.biallelic.het.vcf

for CHR in $(seq 1 22; echo X);
do
    VCFS+=( ${TEMP_WORK}/gatk/${ID}.chr${CHR}.knownSNVs.biallelic.het.vcf )
done

bcftools concat -o $OUT "${VCFS[@]}"

gatk IndexFeatureFile -I $OUT
