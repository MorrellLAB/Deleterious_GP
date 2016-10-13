#!/usr/bin/env python
"""Calculate the derived allele frequency from an allele counts report from
PLINK and the SNP effects table."""

import sys

effects = sys.argv[1]
report = sys.argv[2]

ancestral = {}
with open(effects, 'r') as f:
    for index, line in enumerate(f):
        if index == 0:
            continue
        else:
            tmp = line.strip().split()
            if tmp[0] == '-' or tmp[10] == 'N':
                continue
            else:
                ancestral[tmp[0]] = tmp[10]

#   We will assume that The C1, C2, and C3 reports are the same lengths
dafs = []

with open(report, 'r') as f:
    for index, line in enumerate(f):
        if index == 0:
            continue
        else:
            tmp = line.strip().split()
            snpid = tmp[1]
            if snpid not in ancestral:
                continue
            else:
                a1 = tmp[2]
                a2 = tmp[3]
                anc = ancestral[snpid]
                if anc == a1:
                    der = a2
                    d_col = 5
                elif anc == a2:
                    der = a1
                    d_col = 4
                else:
                    dafs.append('NA')
                    continue
                if int(tmp[4]) + int(tmp[5]) == 0:
                    dafs.append('NA')
                else:
                    daf = float(tmp[d_col]) / (int(tmp[4]) + int(tmp[5]))
                dafs.append(str(daf))

print '\n'.join(dafs)
