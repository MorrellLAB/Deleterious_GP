#!/bin/bash

#PBS -l mem=28gb,nodes=1:ppn=7,walltime=48:00:00
#PBS -m abe
#PBS -M konox006@umn.edu
#PBS -q mesabi

set -e
set -u
set -o pipefail

module load parallel
module load freebayes
module load ogap
module load samtools
module load bamtools

#   Sample params
SAMPLE_LIST="/panfs/roc/scratch/tkono/Genomic_Prediction/Sample_List.txt"
REF_GEN="/panfs/roc/groups/9/morrellp/shared/References/Reference_Sequences/Barley/Morex/barley_RefSeq_v1.0/150831_barley_pseudomolecules.fasta"
OUTDIR="/panfs/roc/scratch/tkono/Genomic_Prediction/Variants"
REGIONS="/panfs/roc/scratch/tkono/Genomic_Prediction/Regions_List.txt"
PROJECT="Genomic_Prediction"

#   Other params
THETA=0.008 # Set Theta
PLOIDY=1 # Set ploidy level
MIN_MAPQ=20 # Set minimum mapping quality
MIN_BASQ=20 # Set minimum base quality
MIN_ALTC=20 # Set minimum alternate count
MIN_ALTF=0.3 # Set minimum alternate fraction


#   Check dependencies
declare -a FREE_DEPENDENCIES=(bamtools ogap bamleftalign samtools freebayes)
for prog in "${FREE_DEPENDENCIES[@]}"; do if ! $(command -v "${prog}" > /dev/null 2> /dev/null); then echo "Failed to find ${prog}, exiting..." >&2; exit 1; fi; done

#   Check samples
if ! [[ -f "${SAMPLE_LIST}" ]]; then echo "Failed to find ${SAMPLE_LIST}, exiting..." >&2; exit 1; fi
for sample in $(<${SAMPLE_LIST}); do if ! [[ -f "${sample}" ]]; then echo "Cannot find sample ${sample}, exiting..." >&2; exit 1; fi; done
if ! [[ -f "${REGIONS}" ]]; then echo "Failed to find regions file ${REGIONS}, exiting..." >&2; exit 1; fi

#   A function that can be parallelized
function CallsByRegion() {
    local samples="${1}"
    local reference="${2}"
    local region="${3}"
    local theta="${4}"
    local ploidy="${5}"
    local mapq="${6}"
    local basq="${7}"
    local altc="${8}"
    local altf="${9}"
    local outdir="${10}"
    local project="${11}"
    local date=$(date +%Y-%m-%d)
    local output="${outdir}/${project}_${region}_${date}.vcf"
    time bamtools merge -region ${region} -list "${samples}" \
     | time ogap -z -R 25 -C 20 -Q 20 -S 0 -f "${reference}" \
     | time bamleftalign -f "${reference}" \
     | time samtools calmd -EAru - "${reference}" 2> /dev/null \
     | time freebayes \
        --fasta-reference "${reference}" \
        --bam-list "${samples}" \
        --region ${region} \
        --theta "${theta}" \
        --ploidy "${ploidy}" \
        --use-reference-allele \
        --min-mapping-quality "${mapq}" \
        --min-base-quality "${basq}" \
        --min-alternate-count "${altc}" \
        --min-alternate-fraction "${altf}" \
        --hwe-priors-off \
        --no-mnps \
        --no-complex \
        --vcf ${output}
}

#   Export the function
export -f CallsByRegion

parallel --verbose CallsByRegion ${SAMPLE_LIST} ${REF_GEN} {} ${THETA} ${PLOIDY} ${MIN_MAPQ} ${MIN_BASQ} ${MIN_ALTC} ${MIN_ALTF} ${OUTDIR} ${PROJECT} :::: ${REGIONS}
