#!/bin/bash -l
#PBS -l "nodes=1:ppn=4,mem=12gb,walltime=6:00:00"
#PBS -m abe
#PBS -M konox006@umn.edu
#PBS -A morrellp
#PBS -W group_list=morrellp

# This is a script to make cross-validation datasets for GEMMA GWAS. We will use
# a prdiction approach where we mask 20% of the phenotypic values, fit a model
# for G->P with the other 80%, then predict the masked portion. We will repeat
# this procedure 20 times.

# Load modules
module load plink/1.90b4.1
module load R/3.5.0

# Define paths
#   Work directory
GEMMA_WORK="/panfs/roc/scratch/konox006/GEMMA"
#   GEMMA binary
GEMMA="/home/morrellp/shared/Software/GEMMA-0.98.1/gemma-0.98.1-linux-static"
#   Prefix for imputed genotypes, in PLINK format
SRC_BED="GP_AP"
#   Phenotypic data
PHEN="${GEMMA_WORK}/Adjusted_Phenotypic_Data_800Lines.csv"
#   Script to add phenotypic data
PHENSCRIPT="${GEMMA_WORK}/Add_Phenotype_to_PED.py"
#   Masking script
MASKSCRIPT="${GEMMA_WORK}/Random_Phenotype_Mask.R"

cd "${GEMMA_WORK}"

# First, generate a list of samples to subset the entire pop to just the lines
# with phenotype data
python "${PHENSCRIPT}" "${PHEN}" "${SRC_BED}.fam" 'Yld' \
    | awk '$NF != "-9"' \
    | tee "GP_Yld.fam" \
    | cut -f 1-2 -d ' ' \
    > "pheno_subset.txt"
# Subset with plink
plink \
    --bfile "${SRC_BED}" \
    --allow-extra-chr \
    --keep "pheno_subset.txt" \
    --make-bed \
    --out "${SRC_BED}_Pheno"

# # Copy the phenotyped fam file to the right name. This is because GEMMA will
# # automatically exclude individuals with missing phenotype data from any
# # calculations.
# cp "GP_Yld.fam" "${SRC_BED}_Pheno.fam"

# # Estimate a relatedness matrix of our phenotyped lines. We need the centred
# # relatedness matrix for the breeding value analyses later.
# "${GEMMA}" \
#     -bfile "${SRC_BED}_Pheno" \
#     -gk 1 \
#     -o "relatedness"

# Then, start our loop here.
for i in {00..20}
do
    # Mask the phenotype file
    Rscript "${MASKSCRIPT}" "GP_Yld.fam" > "${SRC_BED}_Pheno.fam"
    # Fit a Bayesian sparse linear mixed model with GEMMA, using a ridge
    # regression method, rather than MCMC
    "${GEMMA}" \
        -bfile "${SRC_BED}_Pheno" \
        -bslmm 2 \
        -o "yld_${i}"
    # Then predict the missing values
    "${GEMMA}" \
        -bfile "${SRC_BED}_Pheno" \
        -predict 1 \
        -epm "${GEMMA_WORK}/output/yld_${i}.param.txt" \
        -emu "${GEMMA_WORK}/output/yld_${i}.log.txt" \
        -o "yld_pred_${i}"
done
