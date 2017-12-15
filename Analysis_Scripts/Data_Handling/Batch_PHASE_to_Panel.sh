#!/bin/bash
# Run the PHASE_to_Panel.py script in batch mode to convert all chromosomes

# Define the paths
GH="${HOME}/Dropbox/GitHub/Deleterious_GP"
SCRIPT="${GH}/Analysis_Scripts/Data_Handling/PHASE_to_Panel.py"
BOPA_PHYS="${GH}/Data/SNP_Positions/384_Physical.vcf"
PHASEDIR="${GH}/Data/Imputation/PHASE/Phase_ParentsOnly_pairs/"
VCFTOOLS=$(which vcftools)

# Define an array for the chromosomes
CHROMS=(
    "chr1H"
    "chr2H"
    "chr3H"
    "chr4H"
    "chr5H"
    "chr6H"
    "chr7H")

# For each chromosome:
for c in ${CHROMS[@]}
do
    python ${SCRIPT} ${PHASEDIR}/parents_${c}_pairs ${BOPA_PHYS} ${c} > ${c}_panel.vcf
done
