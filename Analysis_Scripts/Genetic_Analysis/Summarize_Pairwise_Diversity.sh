#!/bin/bash
# Summarize pairwise diversity for various partitions of SNPs in the genomic
# prediction experiment. This will use the imputed SNPs.

# Set paths
BASE="/Volumes/LaCie/Genomic_Prediction/Summary_Stats/"
VCF="/Volumes/LaCie/Genomic_Prediction/Summary_Stats/AllChr_AlphaPeel_VCF.vcf"
PARTITIONS="/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations"

# Set functional classes and cycles
CYCLES=("parents" "cycle1" "cycle2" "cycle3" "All")
CLASSES=("Noncoding" "Synonymous" "Nonsynonymous" "Deleterious" "All")

for CYC in ${CYCLES[@]}
do
    for FUNC in ${CLASSES[@]}
    do
         if [[ ${CYC} == "All" ]]
            then
             if [[ ${FUNC} == "All" ]]
                then
                vcftools --vcf ${VCF} --site-pi --out GP_${CYC}_${FUNC}
            else
                vcftools --vcf ${VCF} --snps ${PARTITIONS}/GP_${FUNC}.txt --site-pi --out GP_${CYC}_${FUNC}
            fi
        else
             if [[ ${FUNC} == "All" ]]
                then
                vcftools --vcf ${VCF} --keep ${BASE}/${CYC}.list --site-pi --out GP_${CYC}_${FUNC}
            else
                vcftools --vcf ${VCF} --snps ${PARTITIONS}/GP_${FUNC}.txt --keep ${BASE}/${CYC}.list --site-pi --out GP_${CYC}_${FUNC}
            fi
        fi
    done
done
