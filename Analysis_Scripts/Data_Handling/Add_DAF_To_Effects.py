#!/usr/bin/env python
"""Takes the output from a derived allele frequency table and merges it with
the effects table."""

import sys

daf = {}
with open(sys.argv[1], 'r') as f:
    for line in f:
        tmp = line.strip().split()
        snpid = tmp[0]
        anc_base = tmp[1]
        freq = tmp[2]
        daf[snpid] = [anc_base, freq]

with open(sys.argv[2], 'r') as f:
    for index, line in enumerate(f):
        tmp = line.strip().split()
        if index == 0:
            l = tmp[0:10] + ['Anc_Base', 'DAF'] + tmp[10:]
        else:
            if tmp[0] in daf:
                l = tmp[0:10] + daf[tmp[0]] + tmp[10:]
            else:
                l = tmp[0:10] + ['-', '-'] + tmp[10:]
        print '\t'.join(l)
