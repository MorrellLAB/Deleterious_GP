#!/bin/bash
#PBS -l mem=16gb,nodes=1:ppn=16,walltime=24:00:00
#PBS -m abe
#PBS -M konox006@umn.edu
#PBS -q lab
#PBS -A morrellp

module load plink
module load R

cd /panfs/roc/scratch/konox006/Genomic_Prediction/GCTA
CMD=$(sed "${PBS_ARRAYID}q;d" GCTA_cmds.txt)
eval ${CMD}
