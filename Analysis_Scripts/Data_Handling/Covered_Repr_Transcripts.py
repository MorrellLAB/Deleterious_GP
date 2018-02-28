#!/usr/bin/env python
"""Trim down the representative transcripts to only those covered by the 50x
capture alignment. We need this because bedtools can't handle the long
chromosomes of the reference genome for barley. Takes two arguments:
    1) BED of transcripts to keep
    2) GTF of all transcripts"""

import sys

keep = []
with open(sys.argv[1], 'r') as f:
    for line in f:
        keep.append(line.strip().split())

with open(sys.argv[2], 'r') as f:
    for line in f:
        tmp = line.strip().split()
        # Make a search key from chromosome, start, and stop
        key = [tmp[0], tmp[3], tmp[4]]
        if key in keep:
            print line.strip()
        else:
            continue
