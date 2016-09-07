#!/usr/bin/env python
"""Using the pedigree defined by Tyler and the individual IDs for each sample,
assign full family ID and line IDs in the PED files."""

import sys
import math

ped = {}
with open(sys.argv[1], 'r') as f:
    for index, line in enumerate(f):
        if index == 0:
            continue
        else:
            tmp = line.strip().split(',')
            famid = tmp[1].split('-')[0]
            if famid in ped:
                continue
            else:
                mat = tmp[7]
                pat = tmp[8]
                ped[famid] = (mat, pat)

with open(sys.argv[2], 'r') as f:
    for line in f:
        tmp = line.strip().split()
        famid = tmp[0]
        if famid in ped:
            f1 = famid
            f2 = famid + '-' + tmp[1]
            patid = ped[famid][1]
            matid = ped[famid][0]
            sex = '0'
            phen = '-9'
            genotype = tmp[6:]
            print '\t'.join([
                f1,
                f2,
                patid,
                matid,
                sex,
                phen] + genotype)
