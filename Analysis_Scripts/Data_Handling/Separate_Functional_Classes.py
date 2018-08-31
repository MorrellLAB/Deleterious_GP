#!/usr/bin/env python
"""This script will separate functional classes of SNPs in a PLINK MAP file
based on membership in a supplied SNP list. We do it this way because grep is
really slow. Takes two arguments:
    1) MAP file
    2) SNP list (gzipped)
"""

import sys
import gzip

flt = set()
with gzip.open(sys.argv[2], 'rt') as f:
    for line in f:
        flt.add(line.strip())

with open(sys.argv[1], 'r') as f:
    for line in f:
        tmp = line.strip().split()
        snpid = tmp[1].split(';')[0]
        if snpid in flt:
            print(line.strip())
        else:
            continue
