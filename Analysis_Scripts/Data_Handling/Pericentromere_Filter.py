#!/usr/bin/env python
"""Separate the pericentormeric and the euchromatic SNPs from the effects table.
We will use a simple algorithm to separate them. Writes euchromatic SNPs to
stdout and pericentromeric SNPs to stderr. Takes two arguments:
    1) Pericentromeres BED
    2) Effects table
"""

import sys


def main(peri, eff):
    """Main function."""
    p = {}
    with open(peri, 'r') as f:
        for line in f:
            tmp = line.strip().split()
            p[tmp[0]] = (int(tmp[1]), int(tmp[2]))
    with open(eff, 'r') as f:
        for index, line in enumerate(f):
            if index == 0:
                sys.stdout.write(line)
                sys.stderr.write(line)
            else:
                tmp = line.strip().split('\t')
                chrom = tmp[1]
                pos = int(tmp[2])
                if chrom not in p:
                    continue
                else:
                    if pos > p[chrom][0] and pos <= p[chrom][1]:
                        sys.stderr.write(line)
                    else:
                        sys.stdout.write(line)
    return


if len(sys.argv) != 3:
    print """Separate the pericentormeric and the euchromatic SNPs from the effects table.
We will use a simple algorithm to separate them. Writes euchromatic SNPs to
stdout and pericentromeric SNPs to stderr. Takes two arguments:
    1) Pericentromeres BED
    2) Effects table"""
    exit(1)
else:
    main(sys.argv[1], sys.argv[2])
