#!/usr/bin/env python
"""A very specialized script to slice down the derived allele frequency table
with a subset of SNPs that were identified to be interesting by GEMMA LMM tests
and examining haplotypes. Takes three arguments:
    1) 4H derived states
    2) DAF table
    3) deleterious names (gzipped)
"""

import sys
import gzip

dsnp_names = set()
with gzip.open(sys.argv[3], 'rt') as f:
    for line in f:
        dsnp_names.add(line.strip())

yld_don_reg_names = set()
with open(sys.argv[1], 'r') as f:
    for index, line in enumerate(f):
        if index == 0:
            continue
        else:
            yld_don_reg_names.add(line.strip().split()[2])

# Then, slice down the frequency table based on membership in deleterious names
# being in the 4H region of interest
with open(sys.argv[2], 'r') as f:
    for index, line in enumerate(f):
        if index == 0:
            continue
        else:
            snp_name = line.strip().split()[0]
            if snp_name in dsnp_names and snp_name in yld_don_reg_names:
                print(line.strip())
