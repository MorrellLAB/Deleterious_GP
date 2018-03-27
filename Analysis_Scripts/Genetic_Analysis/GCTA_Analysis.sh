#!/bin/bash
# Script to perform the analysis for GCTA. This is a wrapper for the wrapper,
# because we want to resample at multiple intensities. And multiple traits. This
# is a big ugly script.

# Set paths
#   Base for output
YIELD_BASE="/Volumes/LaCie/Genomic_Prediction/GCTA/Analysis/Yield_GCTA"
DON_BASE="/Volumes/LaCie/Genomic_Prediction/GCTA/Analysis/DON_GCTA"
HEIGHT_BASE="/Volumes/LaCie/Genomic_Prediction/GCTA/Analysis/Height_GCTA"
#   The GCTA control script
GCTA_SCRIPT="/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Analysis_Scripts/Genetic_Analysis/GCTA_Resampling.sh"
#   The frequency distributions for the partitions of SNPs
PARTITIONS=(
    "/Volumes/LaCie/Genomic_Prediction/GCTA/Freqs/noncoding.frq"
    "/Volumes/LaCie/Genomic_Prediction/GCTA/Freqs/synonymous.frq"
    "/Volumes/LaCie/Genomic_Prediction/GCTA/Freqs/nonsynonymous.frq"
    "/Volumes/LaCie/Genomic_Prediction/GCTA/Freqs/deleterious.frq"
)
DELETERIOUS="/Volumes/LaCie/Genomic_Prediction/GCTA/Freqs/deleterious.frq"
#   The basenames of the yield and DON PLINK files
YIELD_PED="/Volumes/LaCie/Genomic_Prediction/GCTA/Source/Yield/Yield"
DON_PED="/Volumes/LaCie/Genomic_Prediction/GCTA/Source/DON/DON"
HEIGHT_PED="/Volumes/LaCie/Genomic_Prediction/GCTA/Source/Height/Height"

# The resampling intensities of the SNPs that we want to do
INTENSITIES=(250 500 750 1000 1250 1500)
# Set the number of replicates to resample for each partition
ITERS="500"

for i in ${INTENSITIES[@]}
do
    for p in ${PARTITIONS[@]}
    do
        # Get the basename of the frequency partition. This will be a directory
        # name
        part=$(basename ${p})
        part_f=${part/.frq/}
        # Start with yield
        #mkdir -p ${YIELD_BASE}/${i}/${part_f}; cd ${YIELD_BASE}/${i}/${part_f}; bash ${GCTA_SCRIPT} ${ITERS} ${YIELD_PED} ${p} ${DELETERIOUS} ${i}
        # Then do DON
        #mkdir -p ${DON_BASE}/${i}/${part_f}; cd ${DON_BASE}/${i}/${part_f}; bash ${GCTA_SCRIPT} ${ITERS} ${DON_PED} ${p} ${DELETERIOUS} ${i}
        # And height
        mkdir -p ${HEIGHT_BASE}/${i}/${part_f}; cd ${HEIGHT_BASE}/${i}/${part_f}; bash ${GCTA_SCRIPT} ${ITERS} ${HEIGHT_PED} ${p} ${DELETERIOUS} ${i}
    done
done
