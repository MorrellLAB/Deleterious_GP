#!/usr/bin/env python
"""A very *very* simple report on whether or not a SNP shows evidence of
segregation distortion. Does so by comparing the allele frequencies in the
parental lines and in their progeny. If SNPs are polymorphic in the parents, but
monomorphic in the progeny, the SNP is said to show segregation distortion.
Uses the frequency outputs from PLINK2."""

import sys
import math

parents = sys.argv[1]
progeny = sys.argv[2]

parental_freqs = {}
progeny_freqs = {}

with open(parents, 'r') as f:
    for index, line in enumerate(f):
        if index == 0:
            continue
        else:
            tmp = line.strip().split()
            if tmp[4] == 'NA':
                continue
            parental_freqs[tmp[1]] = float(tmp[4])

with open(progeny, 'r') as f:
    for index, line in enumerate(f):
        if index == 0:
            continue
        else:
            tmp = line.strip().split()
            if tmp[4] == 'NA':
                continue
            else:
                progeny_freqs[tmp[1]] = float(tmp[4])

for snpid in sorted(parental_freqs.keys()):
    if parental_freqs[snpid] == 0.0:
        flag = 'NA'
    else:
        flag = parental_freqs[snpid] - progeny_freqs[snpid]
    print snpid, flag
