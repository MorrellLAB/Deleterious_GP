#!/bin/bash

#   Script to run PLINK2 (Plink! 1.90) to filter phased and imputed
#   genotyping matrix, run a simple association test, and save the output from
#   the permutation tests. Requires PLINK2, ESMTest and HDF5 utils.
#       PLINK2: https://www.cog-genomics.org/plink2
#       ESMTest: https://github.com/ThorntonLab/ESMtest
#       HDF5 Utils: https://www.hdfgroup.org/downloads/

#   Common PLINK2 options used:
#       --vcf-min-gp: Min GP (genotype probability) to be called not missing
#       --allow-extra-chr: Allow non-numeric (1H, 2H, ...) chromosome names
#       --real-ref-alleles and --keep-allele-order: keep the order of A1 and
#           A2, as defined in the VCF, rather than A1=min A2=maj
#       --allow-no-sex: Don't set ambiguous sex individuals to missing data
#       --nonfounders: Calculate statistics on "founders" and non founders
#           Default is to only work on founders.
#       --r2: Calculate pairwide LD as r-squared
#       --assoc: Run an association mapping analysis
#       --mperm INT: Run INT iterations of the permutation test
#       --merm-save-all: Save all chi-squared stats for each permutation, for
#           each marker

#   Set paths to scripts here
FAM_SCRIPT="${HOME}/DataDisk/Dropbox/Projects/DM_GenomicPrediction/Scripts/Add_Family_to_PED.py"
PHEN_SCRIPT="${HOME}/DataDisk/Dropbox/Projects/DM_GenomicPrediction/Scripts/Add_Phenotype_to_PED.py"

#   Set paths to data files
PEDIGREES="${HOME}/DataDisk/Dropbox/GitHub/Deleterious_GP/Data/Population_Pedigrees.csv"
PHENOTYPES="${HOME}/DataDisk/Dropbox/Projects/DM_GenomicPrediction/Data/Phenotypic_Data/Yield_BLUEs.csv"

#   Which chromosome are we working on?
CHR="$1"
#   Any genotype with a probability less than 0.7 will be set to missing
VCF_MIN_GP="0.7"
#   How many permutations will we do for the association test?
PERMS="1500000"
#   Random seed
SEED="12345"
#   Set the output prefix for the permutations file
PERMOUT="${HOME}/DataDisk/tmp/Perms"
mkdir -p ${PERMOUT}/${CHR}

#   cd to the correct directory
DATADIR="${HOME}/DataDisk/Dropbox/Projects/DM_GenomicPrediction/Data/Progeny_Genotypes/ESM/${CHR}"
cd ${DATADIR}

#   First step: Convert phased VCF to PED/MAP
plink2 \
    --vcf C0C1C2C3_${CHR}.vcf.gz \
    --recode \
    --out plink_${CHR} \
    --vcf-min-gp ${VCF_MIN_GP} \
    --allow-extra-chr \
    --const-fid "0" \
    --real-ref-alleles \
    --keep-allele-order

#   Second step: Add pedigree and phenotype to PED
python ${FAM_SCRIPT} ${PEDIGREES} plink_${CHR}.ped > plink_${CHR}_WithFam.ped
python ${PHEN_SCRIPT} ${PHENOTYPES} plink_${CHR}_WithFam.ped > plink_${CHR}_WithFam_WithPhen.ped

#   Third step: Calculate genotypic frequencies for each marker, and get a list
#   of those that have less than 1 in one of the homozygote classes. We also
#   remove any markers that fall below 1.45% frequency (0.0145, 48 individuals in
#   this study).
plink2 \
    --ped plink_${CHR}_WithFam_WithPhen.ped \
    --map plink_${CHR}.map \
    --freqx \
    --out ${CHR}_freq \
    --allow-extra-chr \
    --allow-no-sex \
    --nonfounders
awk '{ if($5 < 1 || $7 < 1) print $2 }' ${CHR}_freq.frqx > ${CHR}_Excl.txt

plink2 \
    --ped plink_${CHR}_WithFam_WithPhen.ped \
    --map plink_${CHR}.map \
    --freq \
    --out ${CHR}_freq \
    --allow-extra-chr \
    --allow-no-sex \
    --nonfounders
awk '$5 < 0.015 {print $2}' ${CHR}_freq.frq >> ${CHR}_Excl.txt

#   Fourth step: Calculate the missingness on a per-marker basis, and get a
#   list of markers with at least 20% missing calls.
plink2 \
    --ped plink_${CHR}_WithFam_WithPhen.ped \
    --map plink_${CHR}.map \
    --missing \
    --out ${CHR}_miss \
    --allow-extra-chr \
    --allow-no-sex \
    --nonfounders
awk '$5 > 0.20 {print $2}' ${CHR}_miss.lmiss >> ${CHR}_Excl.txt

#   Fifth step: Sort the list of markers to exclude, remove duplicates, and
#   trim the genotyping matrix. Make the binary PLINK files for input into
#   'perms2h5' later.
sort ${CHR}_Excl.txt | uniq > t.txt
mv t.txt ${CHR}_Excl.txt
plink2 \
    --ped plink_${CHR}_WithFam_WithPhen.ped \
    --map plink_${CHR}.map \
    --exclude ${CHR}_Excl.txt \
    --recode \
    --out plink_${CHR}_Trimmed \
    --allow-extra-chr \
    --real-ref-alleles \
    --keep-allele-order
plink2 \
    --file plink_${CHR}_Trimmed \
    --make-bed \
    --out plink_${CHR}_Bin \
    --allow-extra-chr \
    --real-ref-alleles \
    --keep-allele-order

#   Sixth step: Calculate pairwise LD as r-squared
plink2 \
    --file plink_${CHR}_Trimmed \
    --r2 \
    --out ${CHR}_LD \
    --allow-extra-chr \
    --allow-no-sex \
    --nonfounders

#   Seventh step: Run the permutation test with association
plink2 \
    --file plink_${CHR}_Trimmed \
    --assoc \
    --mperm ${PERMS} \
    --mperm-save-all \
    --out ${PERMOUT}/${CHR}/${CHR}_Permutations \
    --seed ${SEED} \
    --allow-extra-chr \
    --allow-no-sex

#   Eighth step: Convert the permutations file to HDF5 format, for input into
#   the ESM test software.
perms2h5 \
    -i ${PERMOUT}/${CHR}/${CHR}_Permutations.mperm.dump.all \
    -o ${PERMOUT}/${CHR}/${CHR}_Permutations.h5 \
    -b plink_${CHR}_Bin.bim \
    -n 50 \
    -l ${CHR}_LD.ld

#   Ninth step: Remove the .mperms.dump.all file. Too much space.
rm ${PERMOUT}/${CHR}/${CHR}_Permutations.mperm.dump.all
rm ${PERMOUT}/${CHR}/*.qassoc.mperm
