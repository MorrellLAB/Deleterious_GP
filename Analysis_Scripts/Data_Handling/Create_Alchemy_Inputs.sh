#!/bin/env bash

#    A script to create the input files for Alchemy genotype calling.
#    Using a script called 'Make_Alchemy_Inputs.py' to create input files.
#    Need to iterate over a set of directories

set -e
set -o pipefail

#    Directory containing 
INTENSITIES_DIR=/Users/emilyvonderharr/Dropbox_Morrell_Lab/Genomic_Predictions/alchemy_genotype_calling/Alchemy_Input/
#    Create an array of intensity files 
INTENSITIES=(find ${INTENSITIES_DIR} -name "GS*Custom_BOPAC.txt")

#    Location of BOPA 384 SNP AB genotype to nucleotide translation table
TRANSLATE=/Users/emilyvonderharr/Documents/Github/Deleterious_GP/Data/Genotyping_Data/BOPA_384_AB_genotypes.csv

#    Specify the path to the Python program
MAKE_ALCHEMY_INPUTS=/Users/emilyvonderharr/Documents/Github/Deleterious_GP/Analysis_Scripts/Data_Handling/Make_Alchemy_Inputs.py

#    Iterate over all of the intensity files
for i in "${INTENSITIES[@]}"
    do
        echo "${i}"
        TMP=$(basename "${i}" .txt)
        python $MAKE_ALCHEMY_INPUTS "${i}" "${TRANSLATE}" "${TMP}"
    done