#!/bin/env bash

#    A script to translate a set of raw intensity files from Genome Studio to a standard 'BOPA_C' SNP names.
#    Uses a Python script 'Translate_GS_SNP_Names.py' and iterates over a list of intensity files.

set -e
set -o pipefail

module load python3_ML/3.6.1

#    A SNP translation table list various SNP names
TS_TABLE=/home/morrellp/vonde026/Alchemy/BOPA_Name_Translation.csv
#    A list of intensity files; may have older SNP names
FILE=/home/morrellp/vonde026/Alchemy/2008_Plate_MN_ND_Sample_List.txt
#    A directory to hold intensity files after SNPs are renamed
OUT_DIR=/home/morrellp/vonde026/Alchemy
#    Specify the path to the Python program
TRANSLATE_GS_SNP_NAMES=/home/morrellp/vonde026/Deleterious_GP/Analysis_Scripts/Data_Handling/Translate_GS_SNP_Names.py
#    Read the file list into an array
GSFP=($(cat "${FILE}"))

#    Iterate over the list of intensity files
for i in "${GSFP[@]}"
do
	TMP=$(basename ${i} .txt)
	python ${TRANSLATE_GS_SNP_NAMES} ${TS_TABLE} ${i} > ${OUT_DIR}/${TMP}_BOPAC.txt
done