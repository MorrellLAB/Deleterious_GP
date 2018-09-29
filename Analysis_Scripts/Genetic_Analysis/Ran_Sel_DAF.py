#!/usr/bin/env python
"""Calculate the derived allele frequency for each SNP in the various subsets
of the random and selected lines. Takes two arguments:
    1) Phenotypic data CSV with cycle/ran/sel designation
    2) Derived allele table (gzipped)
"""

import sys
import gzip

# Define a dictionary that has the partitions of the panel These will be a
# bunch of sets for lookup
partitions = {
    'C0_par': set(),
    'C1_ran': set(),
    'C1_sel': set(),
    'C2_ran': set(),
    'C2_sel': set(),
    'C3_ran': set(),
    'C3_sel': set(),
    'all_ran': set(),
    'all_sel': set()
    }

# Make a similarly-structured dict to hold indices
part_idx = {
    'C0_par': [],
    'C1_ran': [],
    'C1_sel': [],
    'C2_ran': [],
    'C2_sel': [],
    'C3_ran': [],
    'C3_sel': [],
    'all_ran': [],
    'all_sel': []
    }
with open(sys.argv[1], 'r') as f:
    for index, line in enumerate(f):
        if index == 0:
            continue
        else:
            tmp = line.strip().split(',')
            lid = tmp[0]
            part = tmp[3]
            partitions[part].add(lid)
            if 'ran' in part:
                partitions['all_ran'].add(lid)
            elif 'sel' in part:
                partitions['all_sel'].add(lid)

# For each SNP, slice up the genotype frequencies by partition
with gzip.open(sys.argv[2], 'rt') as f:
    for index, line in enumerate(f):
        if index == 0:
            # Split up the line IDs by partition
            header = line.strip().split()
            for p in partitions:
                for i, lid in enumerate(header):
                    if lid in partitions[p]:
                        part_idx[p].append(i)
            # Print the header line
            h = ['Chromome', 'Pos', 'SNP_ID', 'Anc_State']
            for p in sorted(partitions):
                h.append(p)
            print('\t'.join(h))
        else:
            tmp = line.strip().split()
            # If the ancestral state is 'N', we can't calcualte DAFs
            if tmp[3] == 'N' or tmp[4] == 'N':
                dafs = [tmp[0], tmp[1], tmp[2], 'N', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA']
            else:
                dafs = [tmp[0], tmp[1], tmp[2], tmp[3]]
                for part in sorted(partitions):
                    p_der = [tmp[i] for i in part_idx[part]]
                    alleles = ''.join(p_der)
                    der_f = str(alleles.count('D')/len(alleles))
                    dafs.append(der_f)
            print('\t'.join(dafs))
