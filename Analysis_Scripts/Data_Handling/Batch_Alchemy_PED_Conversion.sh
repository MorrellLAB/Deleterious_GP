#!/bin/bash
# Simple script to run a batch conversion of an ALCHEMY report to a PLINK PED
# and MAP combination. We do this because one argument will change, and four
# will stay constant.

# Define the argument as something nice
ALCHEMY_IN="${1}"

# Path to the script that converts ALCHEMY output to PLINK
SCRIPT="/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Analysis_Scripts/Data_Handling/Alchemy_to_PLINK_PhysicalPos.py"

# Set paths to files that are required for the Python script, but will not need
# to change from run to run.
BOPA_PHYS="/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Data/BOPA1_BOPA2_physical.vcf"
GEN_MAP="/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Data/Munoz2011_BOPA_ConsensusMap.csv"
NAME_TRANS="/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Data/BOPA_Name_Translation.csv"
AB_TRANS="/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Data/Genotyping_Data/BOPA_384_AB_genotypes.csv"

# Set stdout (PED) and stderr (MAP) file names
PED_OUT="${ALCHEMY_IN/.txt/_ALCHEMY.ped}"
MAP_OUT="${ALCHEMY_IN/.txt/_ALCHEMY.map}"

# Run the program
python ${SCRIPT} \
    ${ALCHEMY_IN} \
    ${BOPA_PHYS} \
    ${GEN_MAP} \
    ${NAME_TRANS} \
    ${AB_TRANS} \
    > ${PED_OUT} \
    2> ${MAP_OUT}
