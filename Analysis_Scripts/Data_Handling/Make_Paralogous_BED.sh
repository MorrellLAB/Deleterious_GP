#!/bin/bash
# Make a BED file that reports regions of high heterozygosity in a resequencing
# panel of barley lines. Requires Python and bedtools.

# Define paths
VCF="/Volumes/LaCie/Genomic_Prediction/Het_Filtering/GATK_Capture_WithID.vcf.gz"
TMPDIR="/Volumes/LaCie/Genomic_Prediction/Het_Filtering"
EXCLUDE_SCRIPT="/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Analysis_Scripts/Data_Handling/VCF_Exclude_BED.py"
HET_SCRIPT="/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Analysis_Scripts/Data_Handling/Per-Site_Heterozygosity.py"

# Define constants
HIGH_HET="0.5"
MERGE_THRESHOLD="100"

# First, calculate per-site heterozygosity, and save the positions that have
# high heterozygosity
python ${HET_SCRIPT} ${VCF} \
    | awk -v het="${HIGH_HET}" '$2>=het {print $1}' \
    > "${TMPDIR}/High_Het_SNPs.txt"

# Then, grep them out of the VCF, and convert the positions into BED intervals,
# then merge them
python ${EXCLUDE_SCRIPT} ${TMPDIR}/High_Het_SNPs.txt ${VCF} \
    | bedtools merge -d ${MERGE_THRESHOLD} -i - \
    > "${TMPDIR}/High_Het_Regions.bed"
