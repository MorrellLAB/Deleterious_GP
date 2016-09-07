#!/usr/bin/env python
"""Identify potential strand flipping (reverse complement) errors in a PLINK
PED file. It is a very simplistic test, and does not identify flip errors in
SNPs that matter. Still, it helps avoid problems when converting to r/QTL
formats for analysis. Assumes the first two lines in the PED file are the
parents of the rest of the individuals."""

import sys


def parse_ped(ped):
    """Read the PED and return two items: the line metadata, and the genotyping
    matrix."""
    meta = []
    geno = []
    with open(ped, 'r') as f:
        for line in f:
            tmp = line.strip().split()
            #   save the metadata
            meta.append(tmp[0:6])
            #   Then build diploid genotypes
            nsnps = len(tmp[6:])/2
            genotypes = []
            for i in range(0, nsnps):
                dip = tmp[6 + 2*i] + tmp[7 + 2*i]
                genotypes.append(dip)
            geno.append(genotypes)
    return (meta, geno)


def switch_error(ped_column):
    """Return True or False for each SNP in a PED file for whether it is a
    potential switch error or not. True is if the parents are both homozygous
    for the same allele, and the progeny are monomorphic for a different allele.
    """
    #   The first two elements are the parental calls
    parents = ped_column[0:2]
    #   If the parents are different, then we can't do anything, return False
    if len(set(parents)) != 1:
        return False
    else:
        #   Else, the parents are monomorphic for this marker. We then look at
        #   the progeny.
        progeny = ped_column[2:]
        #   If the progeny are monomorphic, then we compare the genotype to the
        #   parents.
        if len(set(progeny)) == 1:
            if progeny[0] in parents:
                #   If the progeny have identical calls to the parents, it all
                #   checks out
                return False
            else:
                return True
        else:
            #   If the parents are monomorphic but the progeny aren't, we have
            #   bigger problems and have to look at the data closer. We know
            #   that it should be filled with parental calls, though.
            return True


def new_geno_matrix(flip, geno):
    """Flip the genotpyes that need to be flipped according to the flags in the
    flip list. Returns a new genotyping matrix with corrected calls."""
    new_geno = []
    for flag, calls in zip(flip, geno):
        if flag:
            new_geno.append((calls[0],) * len(calls))
        else:
            new_geno.append(calls)
    return new_geno


def main(ped):
    """Main function"""
    metadata, geno_matrix = parse_ped(ped)
    #   We have to transpose the genotyping matrix to iterate over markers
    #   rather than individuals
    geno_matrix_t = zip(*geno_matrix)
    to_flip = [switch_error(col) for col in geno_matrix_t]
    #   Iterate through the transposed genotyping matrix and the list of boolean
    #   values, and flip those that need to be flipped.
    flipped = new_geno_matrix(to_flip, geno_matrix_t)
    #   Then transpose it back
    new_geno = zip(*flipped)
    #   Then, print the new PED file
    for m, g in zip(metadata, new_geno):
        #   Separate the genotypes back into separate fields
        dip = [allele for call in g for allele in call]
        toprint = list(m) + list(dip)
        print '\t'.join(toprint)
    return


main(sys.argv[1])
