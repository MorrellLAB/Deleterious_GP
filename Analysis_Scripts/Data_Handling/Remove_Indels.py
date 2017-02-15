#!/usr/bin/env python
"""Super simple script to filter indels/length polymorphisms from a VCF."""

import sys

with open(sys.argv[1], 'r') as f:
    for line in f:
        if line.startswith('#'):
            print line.strip()
        else:
            tmp = line.strip().split('\t')
            ref = tmp[3]
            alt = tmp[4]
            if len(ref) != 1 or len(alt) != 1:
                continue
            else:
                print line.strip()
