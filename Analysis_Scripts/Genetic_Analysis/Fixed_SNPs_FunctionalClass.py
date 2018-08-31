#!/usr/bin/env python
"""Dumb script to count up the SNPs that are fixed (DAF == 1) and lost 
(DAF == 0) among functional classes in cycle 3. Does not take arguments; the
paths are hardcoded."""

import gzip
import sys

C3_FREQS = '/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/Genotype_Freqs/C3_Freq.frqx.gz'
ANCESTRAL = '/Users/tomkono/Dropbox/Large_Project_Files/Genomic_Prediction/Ancestral_State/GP_Ancestral.txt.gz'
NONCODING = '/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Noncoding.names'
SYNONYMOUS = '/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Synonymous.names'
NONSYNONYMOUS = '/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Nonsynonymous.names'
DELETERIOUS = '/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Deleterious.names'

# Parse the lists and store SNP IDs
nonc = set()
with open(NONCODING, 'r') as f:
    for line in f:
        nonc.add(line.strip())

syn = set()
with open(SYNONYMOUS, 'r') as f:
    for line in f:
        syn.add(line.strip())

deleterious = set()
with open(DELETERIOUS, 'r') as f:
    for line in f:
        deleterious.add(line.strip())


nonsyn = set()
with open(NONSYNONYMOUS, 'r') as f:
    for line in f:
        if line.strip() in deleterious:
            continue
        else:
            nonsyn.add(line.strip())

anc_alleles = {}
with gzip.open(ANCESTRAL, 'rt') as f:
    for index, line in enumerate(f):
        if index == 0:
            continue
        else:
            tmp = line.strip().split()
            snpid = tmp[2]
            ancestral = tmp[3]
            anc_alleles[snpid] = ancestral


# Make a matrix to hold the counts. The rows will be:
#   0 - noncoding
#   1 - synonymous
#   2 - nonsynonymous
#   3 - deleterious
# The columns will be
#   0 - lost
#   1 - segregating
#   2 - fixed
counts = [
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0]
    ]

with gzip.open(C3_FREQS, 'rt') as f:
    for index, line in enumerate(f):
        if index == 0:
            continue
        else:
            tmp = line.strip().split()
            snpid = tmp[1]
            if snpid in nonc:
                row = 0
            elif snpid in syn:
                row = 1
            elif snpid in nonsyn:
                row = 2
            elif snpid in deleterious:
                row = 3
            a1 = tmp[2]
            a2 = tmp[3]
            hom_a1 = tmp[4]
            het = tmp[5]
            hom_a2 = tmp[6]
            astate = anc_alleles.get(snpid, 'N')
            if astate == 'N':
                continue
            elif astate == a1:
                if hom_a1 == '0' and het == '0':
                    # derived allele is fixed
                    col = 2
                elif hom_a2 == '0' and het == '0':
                    # ancestral allele is fixed
                    col = 0
                else:
                    # Still segregating
                    col = 1
            elif astate == a2:
                if hom_a1 == '0' and het == '0':
                    col = 0
                elif hom_a2 == '0' and het == '0':
                    col = 2
                else:
                    col = 1
            else:
                continue
            # Append the right number
            counts[row][col] += 1

# Then print out the counts
classes = ['Noncoding', 'Synonymous', 'Nonsynonymous', 'Deleterious']
print('Class Lost Seg Fixed')
for index, row in enumerate(counts):
    print(classes[index], row[0], row[1], row[2])
