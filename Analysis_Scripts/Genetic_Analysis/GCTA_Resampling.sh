#!/bin/bash

#   Script to perform resampling of the markers in the genomic prediction data
#   to estimate mean and variance of proportion of phenotypic data explained
#   by various subsets of markers.

#   Parameters for the resampling
N_ITERS="5000"
#   Note that some of these executable names will be different on other
#   platforms.
PLINK=$(which plink1.9)
GCTA="/home/tom/DataDisk/Software/GCTA/gcta64"
SUBSAMP="/home/tom/DataDisk/Dropbox/GitHub/Deleterious_GP/Analysis_Scripts/Genetic_Analysis/Downsample_SNPs.R"
GCTA_SOURCE="/home/tom/DataDisk/Dropbox/Projects/DM_GenomicPrediction/Data/Progeny_Genotypes/GCTA/Yield/Source_Files/GCTA_Source_WithPheno"

#   The target frequencies file is the first argument
TARGFREQ="$1"
#   The deleterious frequencies file is the second argument
DELFREQ="$2"

#   Uncomment this when doing a real GCTA analysis
cut -f 1,2,6 -d ' ' "${GCTA_SOURCE}".fam > phenotypes.phen

touch 'VG_VP.txt'
while [ $(wc -l < VG_VP.txt) -lt ${N_ITERS} ]; do
    #   Scramble the phenotpyes with respect to the genotypes. Comment this out
    #   when not doing null expectation generation.
    #cut -f 1,2 -d ' ' "${GCTA_SOURCE}".fam > p.txt
    #cut -f 6 -d ' ' "${GCTA_SOURCE}".fam | shuf - > t.txt
    #paste p.txt t.txt > phenotypes.phen
    #   Sample a subset of the target markers according to the frequency
    #   distribution of the deleterious markers
    Rscript "${SUBSAMP}" "${DELFREQ}" "${TARGFREQ}" > chosen.txt
    #   Trim down the large marker list with PLINK
    ${PLINK}\
        --bfile "${GCTA_SOURCE}"\
        --extract chosen.txt\
        --out tmp\
        --allow-no-sex\
        --make-bed
    #   Then, estimate the GRM with the chosen markers
    ${GCTA} --bfile tmp --make-grm-bin --out tmp_GRM
    #   And then estimate the proportion of phenotypic variance explained by
    #   the markers. Set sevral VG/VP priors to make sure we converge on a
    #   good estimate.
    ${GCTA}\
    --reml\
    --reml-no-constrain\
    --grm-bin tmp_GRM\
    --pheno phenotypes.phen\
    --out tmp_VAR\
    --reml-priors 0.05 0.15 0.25 0.35 0.45\
    --reml-maxit 2000
    #   Then get the variance components out of the output
    grep 'V(G)\s' tmp_VAR.hsq | cut -f 2,3 >> VG.txt
    grep 'V(e)' tmp_VAR.hsq | cut -f 2,3 >> VE.txt
    grep '^Vp' tmp_VAR.hsq | cut -f 2,3 >> VP.txt
    grep 'V(G)/Vp' tmp_VAR.hsq | cut -f 2,3 >> VG_VP.txt
    #   Then, clean up!
    rm chosen.txt tmp*
done
