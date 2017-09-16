#!/bin/bash
# Convert a batch of ALCHEMY call files to PLINK PED/MAP files. Requires the
# GNU coreutils on Mac OSX

if [[ "$(uname)" == "Darwin" ]]
then
    CUT=$(which gcut)
    SORT=$(which gsort)
else
    CUT=$(which cut)
    SORT=$(which sort)
fi

# Define paths
GH_ROOT="/Users/tomkono/Dropbox/GitHub/Deleterious_GP"
CONV_SCRIPT="${GH_ROOT}/Analysis_Scripts/Data_Handling/Alchemy_to_PLINK_PhysicalPos.py"
ORDER_SCRIPT="${GH_ROOT}/Analysis_Scripts/Data_Handling/Fix_Map_Order.py"
INTERPOLATE_SCRIPT="${GH_ROOT}/Analysis_Scripts/Data_Handling/Interpolate_MAP_Positions.py"
FAM_SCRIPT="${GH_ROOT}/Analysis_Scripts/Data_Handling/Add_Family_to_PED.py"
PARENT_NAMES="${GH_ROOT}/Analysis_Scripts/Data_Handling/Parental_Names.txt"
SED_PARENTS="${GH_ROOT}/Analysis_Scripts/Data_Handling/Replace_Parent_Names.sed"
ALCHEMY_REPORTS="${GH_ROOT}/Data/Genotyping_Data/Alchemy_Calls/Complete_ALCHEMY"
PHYS_MAP="${GH_ROOT}/Data/BOPA1_BOPA2_physical.vcf"
GEN_MAP="${GH_ROOT}/Data/Munoz2011_BOPA_ConsensusMap.csv"
NAME_TRANS="${GH_ROOT}/Data/BOPA_Name_Translation.csv"
AB_TRANS="${GH_ROOT}/Data/Genotyping_Data/BOPA_384_AB_genotypes.csv"
PEDIGREES="${GH_ROOT}/Data/Population_Pedigrees.csv"

cd ${ALCHEMY_REPORTS}
# Find all files with names that end in 'calls.txt'
REPS=($(find . -type f -name '*calls.txt' | ${SORT} -V))

# Start converting them
for report in ${REPS[@]}
do
    echo "Converting ${report} to PLINK"
    ped_out=${report/_calls.txt/_ALCHEMY.ped}
    map_out=${report/_calls.txt/_ALCHEMY.map}
    python ${CONV_SCRIPT} \
        ${report} \
        ${PHYS_MAP} \
        ${GEN_MAP} \
        ${NAME_TRANS} \
        ${AB_TRANS} \
        > ${ped_out} \
        2> ${map_out}
done

# Next, combine the BOPA1 and BOPA2 SNPs into one file.
BOPA1=($(find . -type f -name 'GS*BOPA1*.ped' | ${SORT} -V))
BOPA2=($(find . -type f -name 'GS*BOPA2*.ped' | ${SORT} -V))
# Iterate through pairs of BOPA1 and BOPA2 files
for i in $(seq 0 $(expr ${#BOPA1[@]} - 1))
do
    # Get the BOPA1 and BOPA2 PED filenames
    bopa1_ped=${BOPA1[$i]}
    bopa2_ped=${BOPA2[$i]}
    # Generate the matching MAP filenames
    bopa1_map=${bopa1_ped/ped/map}
    bopa2_map=${bopa2_ped/ped/map}
    # Extract the year and bp info from the filename
    year=$(echo ${bopa1_ped} | cut -d '_' -f 2 | sed -E 's/BOPA.//g')
    program=$(echo ${bopa1_ped} | cut -d '_' -f 3)
    # Define a combined output filename
    comb_ped="${year}_${program}_BOPA.ped"
    comb_map="${year}_${program}_BOPA.map"
    # Define a sorted output filename
    sorted_map=${comb_map/.map/_Gen.map}
    # Paste the BOPA2 SNPs on to the end of the BOPA1 SNPs. Recall that PED
    # files have 6 columns that we need to exclude from the pasted file. The
    # MAP files can just be concatenated
    echo "Combining ${bopa1_ped} and ${bopa2_ped}"
    ${CUT} -f 1-6 --complement ${bopa2_ped} | paste ${bopa1_ped} - > ${comb_ped}
    cat ${bopa1_map} ${bopa2_map} > ${comb_map}
    # Then, sort the map on genetic position, and reorder the PED to match
    echo "Sorting and interpolating ${comb_map}"
    ${SORT} -k 1,1V -k 3,3V ${comb_map} > ${sorted_map}
    # And interpolate the missing positions
    python ${INTERPOLATE_SCRIPT} ${sorted_map} > GP_Parents.map
    # Reorder the PED, tack it on to a combined PED
    python ${ORDER_SCRIPT} ${comb_map} GP_Parents.map ${comb_ped} >> All_Chips_BOPA.ped
done

# Now, trim down the combined PED to just the genomic prediction population
# parents, and replace their names
echo "Trimming PED to just genomic prediction population parents"
grep -f ${PARENT_NAMES} All_Chips_BOPA.ped | sed -f ${SED_PARENTS} > GP_Parents.ped

# Next, reorder the progeny so that they match the parents
echo "Reordering and concatenating progeny"
python ${ORDER_SCRIPT} Cycle_1_ALCHEMY.map GP_Parents.map Cycle_1_ALCHEMY.ped >> GP_All_Progeny.ped
python ${ORDER_SCRIPT} Cycle_2_ALCHEMY.map GP_Parents.map Cycle_2_ALCHEMY.ped >> GP_All_Progeny.ped
python ${ORDER_SCRIPT} Cycle_3_ALCHEMY.map GP_Parents.map Cycle_3_ALCHEMY.ped >> GP_All_Progeny.ped

# Add family info to the progeny, and exclude lines without pedigree info
echo "Adding pedigree to progeny"
python ${FAM_SCRIPT} ${PEDIGREES} GP_All_Progeny.ped | grep -vE '^0' > GP_Progeny.ped

# Then, combine the parents and progeny
echo "Concatenating parents and progeny"
cat GP_Parents.ped GP_Progeny.ped > GP_All.ped
