#!/bin/env bash

#    A script to translate a set of raw intensity files from Genome Studio to a standard 'BOPA_C' SNP names.
#    Uses a Python script 'Translate_GS_SNP_Names.py' and iterates over a list of intensity files.
 
set -e
set -o pipefail

#    A SNP translation table list various SNP names
TS_TABLE=/Users/emilyvonderharr/Documents/Github/Deleterious_GP/Data/BOPA_Name_Translation.csv

#    A list of intensity files; may have older SNP names
FILE=/Users/emilyvonderharr/Dropbox_Morrell_Lab/Genomic_Predictions/alchemy_genotype_calling/Genome_Studio_Files.txt
#    Read the file list into an array

#    A directory to hold intensity files after SNPs are renamed
OUT_DIR=/Users/emilyvonderharr/Dropbox_Morrell_Lab/Genomic_Predictions/alchemy_genotype_calling/Alchemy_Input

#    Specify the path to the Python program
TRANSLATE_GS_SNP_NAMES=/Users/emilyvonderharr/Documents/Github/Deleterious_GP/Analysis_Scripts/Data_Handling/Translate_GS_SNP_Names.py

#    Iterate over the list of intensity files
for i in "${GSFP[@]}"
do
	TMP=$(basename ${i} .txt)
	python ${TRANSLATE_GS_SNP_NAMES} ${TS_TABLE} ${i} > ${OUT_DIR}/${TMP}_BOPAC.txt
done