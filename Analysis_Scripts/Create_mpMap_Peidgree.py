#!/usr/bin/env python
"""Take the peidgree that Tyler Tiede created, and convert it into a format
that is usable with EH's mpMap package."""

import sys

#   The order of columns in Tyler's pedigree file is:
#       1) index
#       2) Line ID
#       3) Cycle
#       4) Selected/Random
#       5) BP1 x BP2
#       6) ?
#       7) Pedigree string
#       8) Parent 1
#       9) Parent 2

#   The columns for the mpMap pedigree should be
#       1) Line ID ('id')
#       2) Parent 2 ('Male')
#       3) Parent 1 ('Female')
#       4) whether or not the line was observed ('Obs')
#           Not sure what this means, exactly, but I'll fit it in with 1

with open(sys.argv[1], 'r') as f:
    for index, line in enumerate(f):
        if index == 0:
            print '\t'.join(['id', 'Male', 'Female', 'Observed'])
        else:
            tmp = line.strip().split(',')
            print '\t'.join([tmp[1], tmp[9], tmp[8], '1'])
