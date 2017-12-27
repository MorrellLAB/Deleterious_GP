#!/usr/bin/env python
"""Drop in the BOPA SNP IDs for SNPs that have identical position to the exome
capture resequencing variants. Does not replace genotypes. Prints the modified
VCF to stdout, and a log table to stderr. Takes two arguments:
    1) Exome capture VCF (gzipped)
    2) 384 BOPA VCF
"""

import sys
import gzip


def parse_bopa(b):
    """Parse the BOPA VCF and save the physical position. Because of the weird
    way in which we are comparing positions, we will key on a tuple of
    chromosome and physical position."""
    pos = {}
    with open(b, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                tmp = line.strip().split()
                pos[(tmp[0], tmp[1])] = tmp[2]
    return pos


def main(excap, bopa):
    """Main function."""
    bopa_pos = parse_bopa(bopa)
    with gzip.open(excap, 'rb') as f:
        for line in f:
            if line.startswith('#'):
                print line.strip()
            elif line.startswith('chrUn'):
                # Skip the unmapped chromosome.
                continue
            else:
                tmp = line.strip().split()
                if (tmp[0], tmp[1]) in bopa_pos:
                    newid = bopa_pos[(tmp[0], tmp[1])]
                    sys.stderr.write(tmp[2] + '\t' + newid + '\n')
                else:
                    newid = tmp[2]
                print '\t'.join([tmp[0], tmp[1], newid] + tmp[3:])
    return


if len(sys.argv) != 3:
    print """Drop in the BOPA SNP IDs for SNPs that have identical position to the exome
capture resequencing variants. Does not replace genotypes. Prints the modified
VCF to stdout, and a log table to stderr. Takes two arguments:
    1) Exome capture VCF (gzipped)
    2) 384 BOPA VCF"""
    exit(1)
else:
    main(sys.argv[1], sys.argv[2])
