#!/usr/bin/env python
"""Find SNPs that need to be flipped, given a VCF and a plink frequency file.
If the alleles do not match, then print the ID to flip. Takes two arugments:
    1) VCF
    2) Plink frequency report"""

import sys


def read_vcf(v):
    """Read the VCF, and store the alleles and SNP ID in a dictionary."""
    alleles = {}
    with open(v, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                tmp = line.strip().split()
                snpid = tmp[2]
                a1 = tmp[3]
                a2 = tmp[4]
                alleles[snpid] = (a1, a2)
    return alleles


def read_freq(freq):
    """Read the plink frequency report, store the alleles and SNP ID in the same
    format as the VCF alleles."""
    alleles = {}
    with open(freq, 'r') as f:
        for index, line in enumerate(f):
            if index == 0:
                continue
            else:
                tmp = line.strip().split()
                snpid = tmp[1]
                a1 = tmp[2]
                a2 = tmp[3]
                alleles[snpid] = (a1, a2)
    return alleles


def main(vcf, freq):
    """Main function."""
    vcf_snps = read_vcf(vcf)
    frq_snps = read_freq(freq)
    # Then compare the two
    for snp in sorted(vcf_snps):
        v = vcf_snps[snp]
        f = frq_snps[snp]
        if f[0] in v and f[1] in v:
            continue
        else:
            print snp
    return


if len(sys.argv) != 3:
    print """Find SNPs that need to be flipped, given a VCF and a plink frequency file.
If the alleles do not match, then print the ID to flip. Takes two arugments:
    1) VCF
    2) Plink frequency report"""
    exit(1)
else:
    main(sys.argv[1], sys.argv[2])
