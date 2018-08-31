#!/usr/bin/env python
"""Read the ancestral state file and the SNP annotations to make PLINK score
files for estimating mutational burden. We will hard-code the paths to the
files because there are several and supplying them as arguments is a bit too
tedious."""

import sys
import gzip

NONC = '/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Noncoding.txt.gz'
SYN = '/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Synonymous.txt.gz'
NONSYN = '/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Nonsynonymous.txt.gz'
DEL = '/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Deleterious.txt.gz'

ANC = '/Users/tomkono/Dropbox/Large_Project_Files/Genomic_Prediction/Ancestral_State/GP_Ancestral.txt.gz'

# Read them into sets for membership tests
noncoding = set()
with gzip.open(NONC, 'rt') as f:
    for line in f:
        noncoding.add(line.strip())

synonymous = set()
with gzip.open(SYN, 'rt') as f:
    for line in f:
        synonymous.add(line.strip())

nonsynonymous = set()
with gzip.open(NONSYN, 'rt') as f:
    for line in f:
        nonsynonymous.add(line.strip())

deleterious = set()
with gzip.open(DEL, 'rt') as f:
    for line in f:
        deleterious.add(line.strip())

# Remove the deleterious from the nonsynonymous, because these two functional
# classes are not mutually exclusive.
nonsynonymous = nonsynonymous - deleterious

# Then open handles to the score files and write the files as SNPs in the
# various functional groups arise.
noncoding_handle = open('Noncoding_Scores.txt', 'w')
synonymous_handle = open('Synonymous_Scores.txt', 'w')
nonsynonymous_handle = open('Nonsynonymous_Scores.txt', 'w')
deleterious_handle = open('Deleterious_Scores.txt', 'w')

# Iterate through the ancestral state file and write the scores
with gzip.open(ANC, 'rt') as f:
    for index, line in enumerate(f):
        if index == 0:
            continue
        else:
            tmp = line.strip().split()
            gp_id = tmp[2].split(';')[0]
            snpid = tmp[2]
            der = tmp[4]
            # Each allele will have a score of 1, because we are just counting
            score = '1'
            if der == 'N':
                continue
            else:
                towrite = ' '.join([snpid, der, score]) + '\n'
                if gp_id in noncoding:
                    noncoding_handle.write(towrite)
                elif gp_id in synonymous:
                    synonymous_handle.write(towrite)
                elif gp_id in nonsynonymous:
                    nonsynonymous_handle.write(towrite)
                elif gp_id in deleterious:
                    deleterious_handle.write(towrite)
                else:
                    continue

# Clean up
noncoding_handle.flush()
synonymous_handle.flush()
nonsynonymous_handle.flush()
deleterious_handle.flush()

noncoding_handle.close()
synonymous_handle.close()
nonsynonymous_handle.close()
deleterious_handle.close()
