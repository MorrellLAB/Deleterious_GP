#!/usr/bin/env python
"""A simple script to print out the per-site heterozygosity in a VCF file."""

import sys
import gzip

with gzip.open(sys.argv[1], 'rb') as f:
    for line in f:
        if line.startswith('#'):
            continue
        else:
            tmp = line.strip().split()
            genotypes = [g.split(':')[0] for g in tmp[9:]]
            # Count the non-missing genotypes
            non_missing = len(genotypes) - genotypes.count('./.')
            # and print the proportion of hets
            p_het = str(float(genotypes.count('0/1'))/non_missing)
            print tmp[2] + '\t' + p_het
