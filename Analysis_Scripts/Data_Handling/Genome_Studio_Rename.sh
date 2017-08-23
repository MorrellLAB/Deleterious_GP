 #!/bin/bash
 
 set -e
 set -o pipefail

TS_TABLE=/Users/emilyvonderharr/Documents/Github/Deleterious_GP/Data/BOPA_Name_Translation.csv
FILE=/Users/emilyvonderharr/Dropbox_Morrell_Lab/Genomic_Predictions/alchemy_genotype_calling/Genome_Studio_Files.txt 
OUT_DIR=/Users/emilyvonderharr/Dropbox_Morrell_Lab/Genomic_Predictions/alchemy_genotype_calling/Alchemy_Input
GSFP=($(cat ${FILE}))

TRANSLATE_GS_SNP_NAMES=/Users/emilyvonderharr/Documents/Github/Deleterious_GP/Analysis_Scripts/Data_Handling/Translate_GS_SNP_Names.py

for i in "${GSFP[@]}"
do
	TMP=$(basename ${i} .txt)
	python ${TRANSLATE_GS_SNP_NAMES} ${TS_TABLE} ${i} > ${OUT_DIR}/${TMP}_BOPAC.txt
done