#!/usr/bin/env python
"""Identify the sites where the parent and progeny PED files need to be flipped
before merging. These will look like quadallelic sites. Takes two arguments:
    1) Merged PED with problematic sites to identify
    2) MAP file
"""

import sys


def identify_flips(alleles):
    """Use the allele counts vector to identify sites to flip. This is going to
    be a simple check - if all are non-zero, then we want to flip them. We will
    exclude missing calls."""
    non_missing = alleles[:-1]
    # Count how many zeroes are in this list - if it is less than 2, we need to
    # flip them.
    if non_missing.count(0) < 2:
        return True
    else:
        return False


def count_alleles(genotypes):
    """Return a vector of allele counts. The order will be [A, T, C, G, N]."""
    # Define a dictionary for easy counting
    counts = {
        'A': 0,
        'T': 0,
        'C': 0,
        'G': 0,
        '0': 0
        }
    for g in genotypes:
        for a in g:
            counts[a] += 1
    return [counts['A'], counts['T'], counts['C'], counts['G'], counts['0']]


def transpose_ped(ped):
    """Transpose a PED file so we can iterate over markers instead of over
    individuals."""
    tped = []
    with open(ped, 'r') as f:
        for line in f:
            tped.append(line.strip().split()[6:])
    tped = zip(*tped)
    return tped


def snp_order(mapfile):
    """Return a list of SNPs in order."""
    ordered = []
    with open(mapfile, 'r') as f:
        for line in f:
            ordered.append(line.strip().split()[1])
    return ordered


def main(pedfile, mapfile):
    """Main function."""
    ordered_snps = snp_order(mapfile)
    t_ped = transpose_ped(pedfile)
    # Iterate through the PED columns and count the alleles
    for index, snp in enumerate(ordered_snps):
        alleles = count_alleles(t_ped[2*index:2*(index+1)])
        flip = identify_flips(alleles)
        if flip:
            print snp
    pass


if len(sys.argv) != 3:
    print """
Identify the sites where the parent and progeny PED files need to be flipped
before merging. These will look like quadallelic sites. Takes two arguments:
    1) Merged PED with problematic sites to identify
    2) MAP file"""
    exit(1)
else:
    main(sys.argv[1], sys.argv[2])
