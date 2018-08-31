#!/usr/bin/env python
"""Add phenotype data to a PED file. This is used for adding phenotype to a
PED that was generated from VCF, as VCF does not store phenotype."""

import sys

# Dictionary for the columns that are the phenotypes in the CSV
phen_cols = {
    'Yld': 6,
    'DON': 7,
    'Pht': 8
    }

#   Read the phenotype data into a dictionary, keyed on line ID
phenotype = {}
with open(sys.argv[1], 'r') as f:
    for index, line in enumerate(f):
        if index == 0:
            continue
        else:
            #   The phenotypic data files are comma-delimited.
            tmp = line.strip().split(',')
            #   We have to do some fiddling with the line names to get them to
            #   match the nmes in the pedigree files
            lineid = tmp[0]
            #   Some cycle 1 names start with G10W instead of MS10S3
            if lineid.startswith('G10W'):
                lineid = lineid.replace('G10W', 'MS10S3')
            elif lineid.startswith('MS11S3'):
                lineid = lineid.replace('MS11S3', 'MS11S2')
            #   Then, the line ID in the pedigree file has a leading 0, and in
            #   the yield data it does not. But only in some cases.
            if lineid.startswith('MS'):
                if '-' in lineid:
                    parts = lineid.split('-')
                    if len(parts[1]) != 3:
                        lineid = parts[0] + '-0' + parts[1]
            #   If the phenotype column is missing, then we skip it and move on
            #   to the next row.
            if tmp[phen_cols[sys.argv[3]]] == 'NA':
                continue
            else:
                phen = float(tmp[phen_cols[sys.argv[3]]])
                phenotype[lineid] = phen

#   Iterate through the PED file, and drop in the phenotype data.
with open(sys.argv[2], 'r') as f:
    for line in f:
        tmp = line.strip().split()
        lineid = tmp[1]
        if lineid not in phenotype:
            phen = '-9'
        else:
            #   take an arithmetic mean of the observations
            phen = str(phenotype[lineid])
        print(' '.join(tmp[0:5] + [phen] + tmp[6:]))
