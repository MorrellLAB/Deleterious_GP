#!/bin/env bash

#PBS -l mem=16gb,nodes=1:ppn=8,walltime=24:00:00
#PBS -m abe
#PBS -M vonde026@umn.edu
#PBS -q lab
#PBS -e /panfs/roc/groups/9/morrellp/vonde026/Error_Files 
#PBS -o /panfs/roc/groups/9/morrellp/vonde026/Error_Files

#    A script to execute Alchemy genotype calling on a set of samples.
#    Alchemy requires 1) an intensity file, 2) a sample sheet, and 3) a SNP map listing positions

set -e
set -o pipefail

#    Load the version of Alchemy install on MSI
module load alchemy_ML

#    Directory containing intensity files
INTENSITY_DIR=${HOME}/Alchemy/Alchemy_Intensities
INPUT_DIR=${HOME}/Alchemy/Alchemy_Inputs
OUT_DIR=${HOME}/Alchemy/Alchemy_Outputs

#    Alchemy intensities files have a name like "GS_2006BOPA1_BA_Plate7_Custom_BOPAC.txt"
#    Saving the name and path to an array
declare -a INTENSITIES=($(find "${INTENSITY_DIR}" -maxdepth 1 -name "GS_*Custom_BOPAC.txt"))
#    Taking the path name from the 1st entry in the array
INTENSITY_PATH=$(dirname "${INTENSITIES[0]}")

#    Iterate over the array of intensity files, and for each sample, cut down the path to the base name
#    The combination of basename and cut return only "GS_2006BOPA1_BA_Plate7"
#    The path is then built back up to reach each sample
for i in "${INTENSITIES[@]}"
	do
		sample=$(basename "${i}" | cut -d "_" -f 1,2,3,4)
	    	alchemy -f "${INTENSITY_PATH}/${sample}_Custom_BOPAC.txt" -s "${INPUT_DIR}/${sample}_Custom_BOPAC_samp_map.txt" -m "${INPUT_DIR}/${sample}_Custom_BOPAC_snp_map.txt" --log="$OUT_DIR/${sample}_log.txt" --illumina > "${OUT_DIR}/${sample}_calls.txt"
    done
