#!/usr/bin/env python
"""Replace the physical map of SNPs listed in a PLINK MAP file with positions
listed in a VCF. The IDs of the SNPs must match. This script will not check
the identites of the chromosome names. Takes two arguments:
    1) MAP file to replace
    2) VCF with desired positions
"""

import sys

def parse_vcf(vcf):
    """Parse the VCF and store the positions in a dictionary."""
    pos = {}
    with open(vcf, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                tmp = line.strip().split('\t')
                pos[tmp[2]] = tmp[1]
    return pos


def main(plinkmap, vcf):
    """Main function."""
    # First, store the VCF data
    vcf_pos = parse_vcf(vcf)
    # Then iterate throguh the map file
    with open(plinkmap, 'r') as f:
        for line in f:
            tmp = line.strip().split()
            snpid = tmp[1]
            if snpid in vcf_pos:
                newpos = vcf_pos[snpid]
            else:
                newpos = '-9'
            print '\t'.join([tmp[0], tmp[1], tmp[2], newpos])
    return


main(sys.argv[1], sys.argv[2])
