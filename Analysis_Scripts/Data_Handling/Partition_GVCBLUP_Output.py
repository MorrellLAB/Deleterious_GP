#!/usr/bin/env python
"""Partition the GVCBLUP output file based on the SNP IDs given. Takes two
arguments:
    1) SNP IDs
    2) GVCBLUP SNP effects file
"""

import sys

keep = []
with open(sys.argv[1], 'r') as f:
    for line in f:
        keep.append(line.strip())

# Cast to set for speedup
keep = set(keep)

with open(sys.argv[2], 'r') as f:
    for index, line in enumerate(f):
        if index == 0:
            print(line.strip())
        else:
            tmp = line.strip().split()
            if tmp[0] in keep:
                print(line.strip())
            else:
                continue
