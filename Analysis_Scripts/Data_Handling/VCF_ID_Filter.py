#!/usr/bin/env python
"""Simple script to subset a VCF given a list of IDs. Prints the selected
IDs to stdout in VCF format. Takes two arguments:
    1) ID list (Gzipped)
    2) VCF (Gzipped)
"""

import sys
import gzip

keep_snps = set()
with gzip.open(sys.argv[1], 'rt') as f:
    for line in f:
        keep_snps.add(line.strip())

with gzip.open(sys.argv[2], 'rt') as f:
    for line in f:
        if line.startswith('#'):
            print(line.strip())
        else:
            tmp = line.strip().split('\t')
            snpid = tmp[2]
            if snpid in keep_snps:
                print(line.strip())
            else:
                continue
