#!/usr/bin/env python
"""Scan through two VCF files (gzipped), and print a VCF with the sites that
overlap, by position. Does not check samples or genotypes. The header for the
VCF will come from the second file. The first file should be the smaller of the
two, for efficiency reasons. The variant IDs from A will be written into the
ID column of AB.vcf."""


import sys
import gzip


def store_vcf(v):
    """Read through the VCF and return a dictionary that links the chromosome
    and position to a SNP ID."""
    a_dat = {}
    with gzip.open(v, 'rb') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                tmp = line.strip().split()
                chrom = tmp[0]
                pos = tmp[1]
                snp_id = tmp[2]
                a_dat[chrom + ':' + pos] = snp_id
    return a_dat


def main(a_vcf, b_vcf):
    """Main function."""
    a_dat = store_vcf(a_vcf)
    # Iterate through B.vcf.gz
    with gzip.open(b_vcf, 'rb') as f:
        for line in f:
            if line.startswith('#'):
                print line.strip()
            else:
                tmp = line.strip().split()
                # Make a key out of chrom:pos, and check for membership in the
                # A.vcf dict
                k = tmp[0] + ':' + tmp[1]
                if k not in a_dat:
                    continue
                else:
                    toprint = '\t'.join([tmp[0], tmp[1], a_dat[k]] + tmp[3:])
                    print toprint
    return


if len(sys.argv) != 3:
    print """Usage:

Intersect_VCF.py [A.vcf.gz] [B.vcf.gz] > AB.vcf

Scan through two VCF files (gzipped), and print a VCF with the sites that
overlap, by position. Does not check samples or genotypes. The header for the
VCF will come from the second file. The first file should be the smaller of the
two, for efficiency reasons. The variant IDs from A will be written into the
ID column of AB.vcf."""
    exit(1)
else:
    main(sys.argv[1], sys.argv[2])
