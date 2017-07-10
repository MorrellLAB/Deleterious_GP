#!/bin/bash
# Impute SNPs with FSFHap. This script will take a PLINK PED/MAP file pair and
# a pedigree file, and impute genotypes. The data are as follows:
#   - 294 BOPA SNPs (those passing QC from a 384 chip) genotyped on parents
#     and three cycles of a genomic prediction experimental population, 3356
#     individuals total
#   - Consensus map positions from Munoz et al. 2011, and physical map
#     positions generated by Li Lei and Paul Hoffman with the RefSeq 1.0 Morex
#     assembly
#   - Pedigree as defined by Kevin Smith and Tyler Tiede for the spring
#     six-row genomic prediction experiment. There were no backcrosses, and
#     progeny were self-fertilized for two generations to generate F3 progeny.
# Note that the PLINK PED/MAP files need to be sorted by position before they
# are loaded into TASSEL. You can do this by running '--recode' on the data set
# with PLINK 1.9, and having only numeric chromosome names.

# Define paths to the input data
GP_PED="/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Data/Genotyping_Data/Alchemy_Calls/TASSEL_Inputs/GP.ped"
GP_MAP="/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Data/Genotyping_Data/Alchemy_Calls/TASSEL_Inputs/GP.map"
PEDIGREE="/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Data/Genotyping_Data/Alchemy_Calls/TASSEL_Inputs/GP_Tassel_Pedigree_Reduced.txt"

# Define paths to the output data
YMD=$(date +%F)
OUT_DIR="/Volumes/LaCie/Genomic_Prediction/Imputation/FSFHap/${YMD}"
LOG="/Volumes/LaCie/Genomic_Prediction/Imputation/FSFHap/FSFHap_${YMD}.log"
OUT_FNAME="GP_Imputed.hmp.txt"

# Define paths to the TASSEL pipeline script
TASSEL="/Users/tomkono/Software/tassel-5-standalone/run_pipeline.pl"

# Define the parameters for the run
CLUSTER="false"     # Use the cluster algorithm. Works well for F4+
WINDOW_LD="true"    # Use window LD algorithm. Works well for inbred parents.
BC="false"          # Use backcross model.
MULTI_BC="false"    # Use multiple backcross model.
MIN_MAF="0.12"      # Minimum minor allele frequency per family
WINDOW_SIZE="45"    # Number of markers to consider for a haplotype
MIN_R="0.15"         # Minimum correlation (LD) for pruning haplotypes
MAX_MISSING="0.9"   # Maximum missing calls per site, per family
NO_HETS="false"     # Remove heterozygous genotypes
MAX_DIFF="0"        # Number of differences allowed to call haplotypes equal
MIN_HAP="4"         # Minimum obs. for haplotype to be considered valid
OVERLAP="30"        # Overlap between windows
FILL_GAP="true"     # Fill gaps with flanking falues
P_HET="0.1"        # Proportion of heterozygous genotypes
MERGE="false"       # Merge parents/progeny/families
OUT_PARENTS="true"  # Replace missing with parents if no recombination
OUT_NUC="true"      # Not sure. Same description as above in manual
OUT_IUPAC="true"    # Use IUPAC ambiguities for heterozygotes in output

# move to the output directory. For some reason, TASSEL cannot create files
# with meaningful names outside the current working directory.
mkdir -p ${OUT_DIR}
cd ${OUT_DIR}

# Run the command
${TASSEL} \
    -plink \
    -ped ${GP_PED} \
    -map ${GP_MAP} \
    -FSFHapImputationPlugin \
    -pedigrees ${PEDIGREE} \
    -logfile ${LOG} \
    -cluster ${CLUSTER} \
    -windowLD ${WINDOW_LD} \
    -bc ${BC} \
    -multbc ${MULTI_BC} \
    -minMaf ${MIN_MAF} \
    -window ${WINDOW_SIZE} \
    -minR ${MIN_R} \
    -maxMissing ${MAX_MISSING} \
    -nohets ${NO_HETS} \
    -maxDiff ${MAX_DIFF} \
    -minHap ${MIN_HAP} \
    -overlap ${OVERLAP} \
    -fillgaps ${FILL_GAPS} \
    -phet ${P_HET} \
    -merge ${MERGE} \
    -outParents ${OUT_PARENTS} \
    -outNuc ${OUT_NUC} \
    -outIUPAC ${OUT_IUPAC} \
    -endplugin \
    -export

# Then, we need to merge them together. We will need to generate a big ugly
# merge command. It has to be done this way...
#   pipeline.pl -fork1 -h in1.txt -fork2 -h in2.txt ... -forkN -h inN.txt ...
#               -combineN -input1 -input2 ... -inputN -mergeGenotypeTables
#               -export out.hmp.txt -runfork1 -runfork2 ... -runforkN
# We will build this command up iteratively.
# First, get an array of the file names
FNAMES=($(find . -type f -name 'imputed_genotypes*' | sort))
# Then start the command
CMD="${TASSEL} "
# For each filename, put the -forkN -h in.txt bit
COUNTER=1
for f in ${FNAMES[@]}
do
    CMD+="-fork${COUNTER} -h ${f} "
    COUNTER=$((${COUNTER} + 1))
done
# Then add the -combineN bit
CMD+="-combine${COUNTER} "
# Then, add the -inputN bits. We have to reset the counter
COUNTER=1
for f in ${FNAMES[@]}
do
    CMD+="-input${COUNTER} "
    COUNTER=$((${COUNTER} + 1))
done
# And add the -mergeGenotypeTables bit, and the export bit
CMD+="-mergeGenotypeTables -export ${OUT_FNAME} "
# And finally, add the -runforkN bits. Reset the counter again
COUNTER=1
for f in ${FNAMES[@]}
do
    CMD+="-runfork${COUNTER} "
    COUNTER=$((${COUNTER} + 1))
done
# And execute the command
${CMD}