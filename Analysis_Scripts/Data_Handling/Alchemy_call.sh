#!/bin/env bash

#    A script to execute Alchemy genotype calling on a set of samples.
#    Alchemy requires 1) an intensity file, 2) a sample sheet, and 3) a SNP map listing positions

set -e
set -o pipefail

#    Load the version of Alchemy install on MSI
module load alchemy_ML

#    Create an array listing each of our directories
#    Cheating a little, we entered a specific set of plates, will have to change when we genotype other samples!
genotyping_plates=(2006BOPA1_BA_Plate7 2006BOPA1_MN_Plate3 2006BOPA2_BA_Plate7 2006BOPA2_MN_Plate3 2007BOPA1_BA_Plate7 2007BOPA1_MN_Plate3 2007BOPA1_N6_Plate9 2007BOPA2_BA_Plate7 2007BOPA2_MN_Plate3 2007BOPA2_N6_Plate9 2008BOPA1_BA_Plate7 2008BOPA2_BA_Plate7)

#    Iterate over the array of directories

for i in "${genotyping_plates[@]}"
    do
        echo "$i"
        alchemy -f "GS_${i}_Custom.txt" -s "${i}_samp_map.txt" -m "${i}_snp_map.txt" --illumina --log="${i}_log.txt" > ./${i}_calls.txt
    done