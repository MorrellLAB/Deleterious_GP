#!/usr/bin/env python
"""Adds family information into a PED file using the Line ID and the pedigree
file from Tyler. This is used for adding family information to PED files that
are generated from VCF, which do not contain family info."""

import sys

#   Parse the pedigree, saving maternal and paternal parentage
ped = {}
with open(sys.argv[1], 'r') as f:
    for index, line in enumerate(f):
        if index == 0:
            continue
        else:
            tmp = line.strip().split(',')
            famid = tmp[0]
            if famid in ped:
                continue
            else:
                mat = tmp[1]
                pat = tmp[2]
                ped[famid] = (mat, pat)

#   Iterate through the PED file, dropping in the correct maternal and paternal
#   IDs.
with open(sys.argv[2], 'r') as f:
    for line in f:
        tmp = line.strip().split()
        famid = tmp[1].split('-')[0]
        if famid in ped:
            f1 = famid
            f2 = tmp[1]
            patid = ped[famid][1]
            matid = ped[famid][0]
        else:
            f1 = '0'
            f2 = tmp[1]
            patid = '0'
            matid = '0'
        sex = '0'
        phen = tmp[5]
        genotype = tmp[6:]
        print('\t'.join([
            f1,
            f2,
            patid,
            matid,
            sex,
            phen] + genotype))
