#!/bin/bash
# Script to perform the analysis for GCTA. This is a wrapper for the wrapper,
# because we want to resample at multiple intensities. And multiple traits. This
# is a big ugly script.

# Set paths
#   Base for output
YIELD_BASE="/panfs/roc/scratch/konox006/Genomic_Prediction/GCTA/Analysis/Yield_GCTA"
DON_BASE="/panfs/roc/scratch/konox006/Genomic_Prediction/GCTA/Analysis/DON_GCTA"
#   The GCTA control script
GCTA_SCRIPT="/panfs/roc/scratch/konox006/Genomic_Prediction/GCTA/GCTA_Resampling.sh"
#   The frequency distributions for the partitions of SNPs
PARTITIONS=(
    "/panfs/roc/scratch/konox006/Genomic_Prediction/GCTA/Freq/noncoding.frq"
    "/panfs/roc/scratch/konox006/Genomic_Prediction/GCTA/Freq/synonymous.frq"
    "/panfs/roc/scratch/konox006/Genomic_Prediction/GCTA/Freq/nonsynonymous.frq"
    "/panfs/roc/scratch/konox006/Genomic_Prediction/GCTA/Freq/deleterious.frq"
)
DELETERIOUS="/panfs/roc/scratch/konox006/Genomic_Prediction/GCTA/Freq/deleterious.frq"
#   The basenames of the yield and DON PLINK files
YIELD_PED="/panfs/roc/scratch/konox006/Genomic_Prediction/GCTA/Source/Yield"
DON_PED="/panfs/roc/scratch/konox006/Genomic_Prediction/GCTA/Source/DON"

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
        echo "mkdir -p ${YIELD_BASE}/${i}/${part_f}; cd ${YIELD_BASE}/${i}/${part_f}; bash ${GCTA_SCRIPT} ${ITERS} ${YIELD_PED} ${p} ${DELETERIOUS} ${i}"
        # Then do DON
        echo "mkdir -p ${DON_BASE}/${i}/${part_f}; cd ${DON_BASE}/${i}/${part_f}; bash ${GCTA_SCRIPT} ${ITERS} ${DON_PED} ${p} ${DELETERIOUS} ${i}"
    done
done
