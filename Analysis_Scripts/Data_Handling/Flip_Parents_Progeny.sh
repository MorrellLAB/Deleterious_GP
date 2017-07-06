#!/bin/bash
# Script to perform the flipping and dropping of SNPs in the genomic prediction
# genotyping dataset to make it all consistent with the forward strand of the
# Morex reference, as of July 2017.
# Requires PLINK2 and Python2.

# Define paths to directories and necessary data.
BASE_DIR="/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Data/Genotyping_Data/Alchemy_Calls/Flipping"
# Parents PED file is generated from merging then trimming BOPA1+BOPA2 genotype
# experiments from T3. It should be in the same marker order as the progeny
# PED (genetic map order)
PARENT_PED="${BASE_DIR}/Parents.ped"
# Progeny PED file is generated from ALCHEMY reports.
PROGENY_PED="${BASE_DIR}/Progeny.ped"
# The full 377-SNP MAP file in genetic map order
GEN_MAP="${BASE_DIR}/Genomic_Prediction_Combined_MapOrder.map"
# SNPs that need to be flipped in the progeny before merging with the parents.
# These are SNPs that have flipped strand between the 384 custom array that was
# used to genotype the progeny and the BOPA1+BOPA2 chips.
PROG_FLIP="${BASE_DIR}/Progeny_Flip_List.txt"
# SNPs that need to be flipped in the whole population so that the alleles are
# on the forward reference strand.
PAR_FLIP="${BASE_DIR}/BOPA_genotype_RC.txt"
# SNPs that need to be dropped because they have inconsistent mapping locations
# or do not genotype reliably.
DROP="${BASE_DIR}/BOPA_genotype_drop.txt"

# Define the names of the executables and scripts
PLINK=$(which plink2)
# Note that we need the -V option, which is only in GNU sort. On OSX, you can
# install this with 'brew install coreutils' and access it with 'gsort'
SORT=$(which gsort)
REORDER_SCRIPT="/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Analysis_Scripts/Data_Handling/Fix_Map_Order.py"

cd ${BASE_DIR}

# The first step is to flip the progeny PED file
${PLINK} \
    --ped ${PROGENY_PED} \
    --map ${GEN_MAP} \
    --allow-extra-chr \
    --flip ${PROG_FLIP} \
    --recode \
    --out Progeny_Flipped

# Then, we need to re-order the flipped progeny PED file because PLINK sortes
# on physical position
python ${REORDER_SCRIPT} Progeny_Flipped.map ${GEN_MAP} Progeny_Flipped.ped > Progeny_Flipped_Reordered.ped

# Concatenate the parents and the flipped-reordered progeny
cat ${PARENT_PED} Progeny_Flipped_Reordered.ped > Combined.ped

# Flip the SNPs that need to be RCed, and drop the ones that are problematic
${PLINK} \
    --ped Combined.ped \
    --map ${GEN_MAP} \
    --allow-extra-chr \
    --flip ${PAR_FLIP} \
    --exclude ${DROP} \
    --recode \
    --out Combined_Flipped

# Sort the genetic map on genetic position
${SORT} -k1,1V -k3,3V Combined_Flipped.map > Combined_Flipped_Gen.map

# Reorder the PED file to match the genetic map order. The
# Combined_Flipped_Gen.ped and Combined_Flipped_Gen.map files are the final
# files from this script. The next steps are to add in the genotyping data
# from M122 and M138, which were not genotyped on the BOPA chips.
python ${REORDER_SCRIPT} Combined_Flipped.map Combined_Flipped_Gen.map Combined_Flipped.ped > Combined_Flipped_Gen.ped
