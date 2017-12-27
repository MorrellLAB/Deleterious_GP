#!/usr/bin/env python
"""This is a really dumb script to make a bunch of 1-bp intervals from a list
of SNP IDs to exclude and a VCF. This is because grep is VERY SLOW for more
than just a couple hundred search strings. Takes two arguments:
    1) SNP IDs to exclude
    2) VCF (gzipped)"""

import sys
import gzip

exclude = []
with open(sys.argv[1], 'r') as f:
    for line in f:
        exclude.append(line.strip())

with gzip.open(sys.argv[2], 'rb') as f:
    for line in f:
        if line.startswith('#'):
            continue
        else:
            tmp = line.strip().split()
            if tmp[2] in exclude:
                chrom = tmp[0]
                start = str(int(tmp[1]) - 1)
                end = tmp[1]
                print chrom + '\t' + start + '\t' + end
            else:
                continue
