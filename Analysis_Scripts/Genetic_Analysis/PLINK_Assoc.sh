#!/bin/bash
# Run a quantitative association analysis with PLINK.

# The name of the PLINK executable
PLINK=$(which plink2)
# The directory where the SNP functional class lists are
SNPNAME_DIR="${HOME}/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations"
# The base directory of the PLINK BED/BIM/FAM with phenotype
PLINK_BASE="${HOME}/Dropbox/Large_Project_Files/Genomic_Prediction/PLINK_qAssoc"

# PLINK filtering criteria
MAX_MAF="0.1"
MIN_MAF="0.02"

for trait in Yield DON Height
do
    for part in Noncoding Synonymous Nonsynonymous Deleterious
    do
        mkdir -p "${PLINK_BASE}/${trait}/${part}"
        plink2 \
            --bfile "${PLINK_BASE}/${trait}/${trait}" \
            --extract "${SNPNAME_DIR}/GP_${part}.names" \
            --allow-extra-chr \
            --allow-no-sex \
            --nonfounders \
            --assoc fisher perm \
            --max-maf "${MAX_MAF}" \
            --maf "${MIN_MAF}" \
            --out "${PLINK_BASE}/${trait}/${part}/${part}_Assoc"
    done
done
