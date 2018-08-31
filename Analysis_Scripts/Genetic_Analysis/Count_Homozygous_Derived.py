#!/usr/bin/env python
"""Count the number of homozygous derived SNPs in each sample in various
functional categories."""

import gzip
import sys

ancestral = '/Users/tomkono/Dropbox/Large_Project_Files/Genomic_Prediction/Ancestral_State/GP_Ancestral.txt.gz'
nc = '/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Noncoding.names'
syn = '/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Synonymous.names'
nonsyn = '/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Nonsynonymous.names'
delet = '/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/SNP_Annotations/GP_Deleterious.names'

noncoding = set()
synonymous = set()
nonsynonymous = set()
deleterious = set()

with open(nc, 'r') as f:
    for line in f:
        noncoding.add(line.strip())

with open(syn, 'r') as f:
    for line in f:
        synonymous.add(line.strip())

with open(nonsyn, 'r') as f:
    for line in f:
        nonsynonymous.add(line.strip())

with open(delet, 'r') as f:
    for line in f:
        deleterious.add(line.strip())

nonsynonymous = nonsynonymous - deleterious

homozygous = {}

with gzip.open(ancestral, 'rt') as f:
    for index, line in enumerate(f):
        if index == 0:
            samples = line.strip().split('\t')[5:]
            for s in samples:
                # Count vector for noncoding, synonymous, nonsynonymous, deleterious
                # respectively.
                homozygous[s] = [0, 0, 0, 0]
        else:
            tmp = line.strip().split()
            snpid = tmp[2]
            anc = tmp[3]
            if anc == 'N':
                continue
            sys.stderr.write(snpid + '\n')
            genotypes = tmp[5:]
            if snpid in noncoding:
                for s, gt in zip(samples, genotypes):
                    if gt == 'DD':
                        homozygous[s][0] += 1
            elif snpid in synonymous:
                for s, gt in zip(samples, genotypes):
                    if gt == 'DD':
                        homozygous[s][1] += 1
            elif snpid in nonsynonymous:
                for s, gt in zip(samples, genotypes):
                    if gt == 'DD':
                        homozygous[s][2] += 1
            elif snpid in deleterious:
                for s, gt in zip(samples, genotypes):
                    if gt == 'DD':
                        homozygous[s][3] += 1
            else:
                continue

# Print it out
print('LineID Noncoding Synonymous Nonsynonymous Deleterious')
for s in sorted(homozygous):
    toprint = ' '.join([
        s,
        str(homozygous[s][0]),
        str(homozygous[s][1]),
        str(homozygous[s][2]),
        str(homozygous[s][3])])
    print(toprint)
