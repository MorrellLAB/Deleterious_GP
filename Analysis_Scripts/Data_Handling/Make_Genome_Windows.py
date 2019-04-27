#!/usr/bin/env python
"""Dumb script to make non-overlapping windows across the barley genome. The
windows size will be given in bp as an argument."""

import sys

try:
    winlen = int(sys.argv[1])
except ValueError:
    sys.stderr.write('Please provide an integer window size.\n')
    exit(2)
except IndexError:
    sys.stderr.write('Dumb script to make non-overlapping windows across the barley genome. The windows size will be given in bp as an argument.\n')
    exit(1)

# Define the chromosome sizes
CHR_SIZES = [
    ('chr1H', 558535432),
    ('chr2H', 768075024),
    ('chr3H', 699711114),
    ('chr4H', 647060158),
    ('chr5H', 670030160),
    ('chr6H', 583380513),
    ('chr7H', 657224000)
    ]

for chrom in CHR_SIZES:
    curr = 0
    while curr < chrom[1]:
        if curr + winlen > chrom[1]:
            end = chrom[1]
        else:
            end = curr + winlen
        print('\t'.join([chrom[0], str(curr), str(end)]))
        curr += winlen
