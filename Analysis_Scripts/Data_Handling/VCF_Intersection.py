#!/usr/bin/env python
"""Dumb script to trim down a VCF given a list of chrom:pos positions. This is
because bedtools cannot intersect large positions, even if it can merge them.
Odd. Takes two arguments:
    1) Positions file
    2) VCF to filter
"""

import sys
import gzip


def main(pos, vcf):
    """Main function."""
    # Read the list of positions to filter
    handle = open(sys.argv[1], 'r')
    pos_flt = [s.strip() for s in handle.readlines()]
    handle.close()
    # Then, iterate through the VCF and print it out
    with gzip.open(vcf, 'rb') as f:
        for line in f:
            if line.startswith('#'):
                print line.strip()
            else:
                p = ':'.join(line.strip().split('\t')[0:2])
                if p in pos_flt:
                    print line.strip()
                else:
                    continue
    return


if len(sys.argv) != 3:
    print """Dumb script to trim down a VCF given a list of chrom:pos positions. This is
because bedtools cannot intersect large positions, even if it can merge them.
Odd. Takes two arguments:
    1) Positions file
    2) VCF to filter"""
else:
    main(sys.argv[1], sys.argv[2])
