#!/usr/bin/env python
"""Simple script to add SNP IDs to a VCF file. Does not read or handle rs IDs.
Takes two arguments:
    1) The VCF to add IDs to
    2) The prefix for the ID. Usually a project name
"""

import sys

vcf = sys.argv[1]
prefix = sys.argv[2]

variants = 1
with open(vcf, 'r') as f:
    for line in f:
        if line.startswith('#'):
            print line.strip()
        else:
            sid = prefix + '_' + str(variants)
            tmp = line.strip().split()
            print tmp[0] + '\t' + tmp[1] + '\t' + sid + '\t' + '\t'.join(tmp[3:])
            variants += 1
