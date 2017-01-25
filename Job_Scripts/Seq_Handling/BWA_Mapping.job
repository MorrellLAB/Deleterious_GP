#!/bin/sh
#PBS -l mem=16gb,nodes=1:ppn=8,walltime=8:00:00 
#PBS -m abe 
#PBS -M konox006@umn.edu 
#PBS -q mesabi

#	Date, for versioning
YMD=$(date +%Y-%m-%d)

#   Programs and versions
#	BWA was cloned on 2016-07-19
BWA=${HOME}/Soft/bwa/bwa
SAMTOOLS=${HOME}/Shared/Software/samtools-1.3.1/samtools
REFERENCE=${HOME}/Shared/References/Reference_Sequences/Barley/Morex/barley_RefSeq_v1.0/barley_pseudomolecules_parts.fa

#   We will use a PBS Task Array to handle batch submission of these alignment
#   jobs. It is a nicer way to submit groups of jobs without having to use
#   a for loop. To populate our array, we will get the contents of the reads
#   directory
READS_DIR=/panfs/roc/scratch/tkono/Genomic_Prediction/Reads
OUTPUT_DIR=/panfs/roc/scratch/tkono/Genomic_Prediction/Aligned
#   Get just the sample names and store them in an array. $() does command
#   substitution, and ($()) will convert that output into an array that we can
#   fetch data from.
#   Sample names are the first field before an underscore.
SAMPLE_NAMES=($(/bin/ls -1 ${READS_DIR} | grep 'Scythe' | cut -f 1 -d '_' | uniq))
#   The sample to analyze is the value of the ${PBS_ARRAYID} variable
SAMPLE=${SAMPLE_NAMES[${PBS_ARRAYID}]}
#   Build the paths to forward and reverse reads from it
FWD=${READS_DIR}/${SAMPLE}_R1_ScytheTrimmed.fastq.gz
REV=${READS_DIR}/${SAMPLE}_R2_ScytheTrimmed.fastq.gz

#	Print some nice metadata for this job
echo "Node list stored in ${PBS_NODEFILE}"
echo "Cores (servers) allocated to this job:"
cat ${PBS_NODEFILE}
echo "Aligning ${SAMPLE}, with task ID of ${PBS_ARRAYID}."
echo "Forward reads in ${FWD}"
echo "Reverse reads in ${REV}"
echo "Output stored in ${OUTPUT_DIR}/${SAMPLE}_${YMD}.bam"

#   Now we run the program with all our options
#   Changes:
#       2016-12-26: Align against the pseudomolecule parts, since it kills GATK
#       2015-12-11: Borrow alignment parameters from Bad Mutations 1 barley
#                   cultivars.
#	Options summary:
#       mem: use the MEM (minimal exact match) algor. exact seed match + SW ext.
#       -t 8: Use 8 threads
#       -k 8: seed length of 8 bases
#       -r 1.0: re-seed of MEM is greater than 1xseed_length
#       -M: mark split hits as secondary alignments, for Picard compat.
#       -T 85: don't print alignments with score less than 85. Match score is
#       1, so the max score for a read is 100. Mismatch penalty is 4, so a 
#       cutoff of 85 filters alignments with more than 3 mismatches in one read.
#       -O 8: Gap open penalty of 8
#       -E 1: Gap extend penalty of 1
#           With no mismatches, these gapping parameters allow an indel of at
#           most 7bp in a single read.
${BWA} mem -t 8 -k 8 -r 1.0 -M -T 85 -O 8 -E 1\
    ${REFERENCE}\
    ${FWD} \
    ${REV} >\
    ${OUTPUT_DIR}/${SAMPLE}_${YMD}.sam
