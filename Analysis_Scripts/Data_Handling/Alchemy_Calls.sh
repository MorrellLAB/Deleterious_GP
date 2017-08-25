#!/bin/env bash

#    A script to execute Alchemy genotype calling on a set of samples.
#    Alchemy requires 1) an intensity file, 2) a sample sheet, and 3) a SNP map listing positions

set -e
set -o pipefail

#    Load the version of Alchemy install on MSI
module load alchemy_ML

#    Directory containing intensity files
INTENSITY_DIR=${HOME}/Alchemy/Alchemy_Intensities/
INPUT_DIR=${HOME}/Alchemy/Alchemy_Inputs
OUT_DIR=${HOME}/Alchemy/Alchemy_Outputs

declare -a INTENSITIES=($(find "${INTENSITY_DIR}" -name "GS*Custom_BOPAC.txt"))
INTENSITY_PATH=$(dirname "${INTENSITIES[0]}")

for i in "${INTENSITIES[@]}"
	do
		sample=$(basename "${i}" | cut -d _ -f 1,2,3,4)
	    	alchemy -f "${INTENSITY_PATH}/${sample}" -s "${INPUT_DIR}/${sample}_samp_map.txt" -m "${INPUT_DIR}/{sample}_snp_map.txt" --log="$OUT_DIR/${sample}_log.txt" --illumina > "${OUT_DIR}/${sample}_calls.txt"
    done
