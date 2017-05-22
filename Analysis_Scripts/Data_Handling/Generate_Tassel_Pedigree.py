#!/usr/bin/env python
"""Convert the PLINK .fam file into a pedigree that works for TASSEL. The
columns for the TASSEL pedigree are as follows:
    Family
    Name
    Parent1
    Parent2
    Contribution1
    Contribution2
    F"""

import sys

#   Keep track of which families we have already printed parental info for
parents = []

with open(sys.argv[1], 'r') as f:
    for index, line in enumerate(f):
        tmp = line.strip().split()
        family = tmp[0]
        if family == '0':
            continue
        lineid = tmp[1]
        p1 = tmp[2]
        p2 = tmp[3]
        contribution = '0.5'
        inbreed = '0.75'
        #   Print out the parents.
        if family not in parents:
            print '\t'.join([
                family,
                tmp[2],
                tmp[2],
                'NA',
                '1',
                '0',
                '1'
                ])
            print '\t'.join([
                family,
                tmp[3],
                tmp[3],
                'NA',
                '1',
                '0',
                '1'
                ])
            parents.append(family)
        print '\t'.join([
            family,
            lineid,
            p1,
            p2,
            contribution,
            contribution,
            inbreed]
            )
