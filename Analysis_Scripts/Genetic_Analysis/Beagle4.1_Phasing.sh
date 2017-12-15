#!/bin/bash
# Use beagle 4.1 to phase and impute the BOPA SNPs in the parental and progeny
# cycles. Note that this requires JRE v8 or higher.

# Define paths
BEAGLE="/Volumes/Scratch/Software/beagle.08Jun17.d8b.jar"

# Directories and required files
GH="/Users/tomkono/Dropbox/GitHub/Deleterious_GP"
PANEL="${GH}/Data/Imputation/Phased_Panel.vcf"
MAP="${GH}/Data/Imputation/Beagle_Inteprolated_MAP.map"
VCF_IN="${GH}/Data/Imputation/Beagle_4.1/BOPA_Beagle_In.vcf"
OUT_DIR="${GH}/Data/Imputation/Beagle_4.1/"

# Define the Beagle options here
OUT_FILE="${OUT_DIR}/Beagle_BOPA_Out"
WINDOW="10"
OVERLAP="4"
SEED="12345"
NITER="10"
IMPUTE="true"
GPROBS="true"
NE="100000"
IBD="false"
IBDTRIM="5"

# Run it
java -jar ${BEAGLE} \
    gt=${VCF_IN} \
    ref=${PANEL} \
    out=${OUT_FILE} \
    map=${MAP} \
    window=${WINDOW} \
    overlap=${OVERLAP} \
    seed=${SEED} \
    niterations=${NITER} \
    impute=${IMPUTE} \
    gprobs=${GPROBS} \
    ne=${NE} \
    ibd=${IBD} \
    ibdtrim=${IBDTRIM}
