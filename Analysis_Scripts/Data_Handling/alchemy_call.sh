#!/bin/env bash

module load alchemy

array=(2006BOPA1_BA_Plate7 2006BOPA1_MN_Plate3 2006BOPA2_BA_Plate7 2006BOPA2_MN_Plate3 2007BOPA1_BA_Plate7 2007BOPA1_MN_Plate3 2007BOPA1_N6_Plate9 2007BOPA2_BA_Plate7 2007BOPA2_MN_Plate3 2007BOPA2_N6_Plate9 2008BOPA1_BA_Plate7 2008BOPA2_BA_Plate7)

for i in "${array[@]}"
    do
        echo "$i"
        alchemy -f "GS_${i}_Custom.txt" -s "${i}_samp_map.txt" -m "${i}_snp_map.txt" --illumina --log="${i}_log.txt" > ./${i}_calls.txt
    done