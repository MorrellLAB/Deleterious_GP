#!/bin/bash
# Automate the PLINK and Python script workflow to combine the parental
# exome capture resequencing VCF with the progeny Alchemy calls.

# Define paths
GH_ROOT="/Users/tomkono/Dropbox/GitHub/Deleterious_GP"
#   The files
PARENTAL_VCF="${GH_ROOT}/Data/SNP_Positions/384-Capture_Overlap.vcf"
ALCHEMY_PED="${GH_ROOT}/Data/Genotyping_Data/Alchemy_Calls/GP_All.ped"
ALCHEMY_MAP="${GH_ROOT}/Data/Genotyping_Data/Alchemy_Calls/GP_All.map"
GEN_MAP="${GH_ROOT}/Data/SNP_Positions/Munoz2011_BOPA_ConsensusMap.csv"
PEDIGREES="${GH_ROOT}/Data/Pedigrees/Population_Pedigrees.csv"
TOFLIP="${GH_ROOT}/Data/Genotyping_Data/BOPA_To_Flip.txt"
#   The scripts
GENMAP_SCRIPT="${GH_ROOT}/Analysis_Scripts/Data_Handling/Add_Genetic_Positions_to_Map.py"
INTERPOLATE_SCRIPT="${GH_ROOT}/Analysis_Scripts/Data_Handling/Interpolate_MAP_Positions.py"
FAM_SCRIPT="${GH_ROOT}/Analysis_Scripts/Data_Handling/Add_Family_to_PED.py"
ORDER_SCRIPT="${GH_ROOT}/Analysis_Scripts/Data_Handling/Fix_Map_Order.py"
FINDFLIP="${GH_ROOT}/Analysis_Scripts/Data_Handling/Find_SNPs_to_Flip.py"

# Convert the parental resequencing to PED/MAP
plink2 \
    --vcf ${PARENTAL_VCF} \
    --allow-extra-chr \
    --keep-allele-order \
    --recode \
    --out Parents_VCF

# Add genetic map information
python \
    ${GENMAP_SCRIPT} \
    Parents_VCF.map \
    ${GEN_MAP} \
    > Parents_VCF_WithGen.map

# Sort the map and interpolate
gsort -k1,1V -k3,3V Parents_VCF_WithGen.map > Parents_VCF_WithGen_Sorted.map
python \
    ${INTERPOLATE_SCRIPT} \
    Parents_VCF_WithGen_Sorted.map \
    > Parents_VCF_WithGen_Sorted_Interpolated.map

# Reorder the PED with the interpolated map
python \
    ${ORDER_SCRIPT} \
    Parents_VCF.map \
    Parents_VCF_WithGen_Sorted_Interpolated.map \
    Parents_VCF.ped \
    > Parents_VCF_Reordered.ped

# Subset the GP_All.ped to only have progeny
awk '$1 != "-9" {print}' ${ALCHEMY_PED} > Progeny.ped

# Subset the progeny to only have excap markers
cut -f 2 Parents_VCF_WithGen_Sorted_Interpolated.map > Keep.txt
plink2 \
    --ped Progeny.ped \
    --map ${ALCHEMY_MAP} \
    --allow-extra-chr \
    --extract Keep.txt \
    --recode \
    --out Progeny_Reseq

# Generate the allele frequency report for the progeny
plink2 \
    --file Progeny_Reseq \
    --allow-extra-chr \
    --freq \
    --nonfounders \
    --out Progeny_Freq

# Flip SNPs
plink2 \
    --file Progeny_Reseq \
    --allow-extra-chr \
    --flip ${TOFLIP} \
    --recode \
    --out Progeny_Flipped

# Reorder flipped progeny
python \
    ${ORDER_SCRIPT} \
    Progeny_Flipped.map \
    Parents_VCF_WithGen_Sorted_Interpolated.map \
    Progeny_Flipped.ped \
    > Progeny_Flipped_Reordered.ped

# Concatenate
cat Parents_VCF_Reordered.ped Progeny_Flipped_Reordered.ped > Combined.ped

# Add family data
python \
    ${FAM_SCRIPT} \
    ${PEDIGREES} \
    Combined.ped \
    > Combined_WithPedigree.ped
