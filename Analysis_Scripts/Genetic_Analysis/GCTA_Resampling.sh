#!/bin/bash

#   Script to perform resampling of the markers in the genomic prediction data
#   to estimate mean and variance of proportion of phenotypic data explained
#   by various subsets of markers.

#   Parameters for the resampling
N_ITERS="5000"
NMARKERS="101"
#   Note that some of these executable names will be different on other
#   platforms.
SHUF=$(which shuf)
SORT=$(which sort)
PLINK=$(which plink1.9)
GCTA="/home/tom/DataDisk/Software/GCTA/gcta64"

#   The PLINK input prefix is the first argument
PLINKFILE="$1"

for ((i=1;i<=${N_ITERS};i++)) do
    #   Randomly sample a subset of the markers
    cut -f 2 ${PLINKFILE}.bim | ${SHUF} | head -n ${NMARKERS} | ${SORT} -V | tee chosen.txt | tr '\n' ' ' >> Sampled_Markers.txt
    #   We need to put a newline into the file, since tr will remove them.
    echo "" >> Sampled_Markers.txt
    #   Trim down the large marker list with PLINK
    ${PLINK} --bfile ${PLINKFILE} --extract chosen.txt --out tmp --allow-no-sex --make-bed
    #   Then, estimate the GRM with the chosen markers
    ${GCTA} --bfile tmp --make-grm-bin --out tmp_GRM
    #   And then estimate the proportion of phenotypic variance explained by
    #   the markers
    ${GCTA} --reml --grm-bin tmp_GRM --pheno ${PLINKFILE}.phen --out tmp_VAR
    #   Then get the variance components out of the output
    grep 'V(G)\s' tmp_VAR.hsq | cut -f 2,3 >> VG.txt
    grep 'V(e)' tmp_VAR.hsq | cut -f 2,3 >> VE.txt
    grep '^Vp' tmp_VAR.hsq | cut -f 2,3 >> VP.txt
    grep 'V(G)/Vp' tmp_VAR.hsq | cut -f 2,3 >> VG_VP.txt
    #   Then, clean up!
    rm chosen.txt tmp*
done
