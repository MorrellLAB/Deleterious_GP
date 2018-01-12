#!/usr/bin/env python
"""Summarize the pre-imputation genotype file that has the BOPA SNPs. We do it
this way because we want to keep track of the number of segregating sites that
are expected within each family. Takes two arguments:
    1) AlphaPeel BOPA input file
    2) Pedigree
"""

import sys


def parse_ped(p):
    """Parse the pedigree. For each individual, store its parents."""
    ped = {}
    with open(p, 'r') as f:
        for line in f:
            tmp = line.strip().split()
            ped[tmp[0]] = (tmp[1], tmp[2])
    return ped


def parse_geno(g):
    """Parse the genotypes and store them for each individual."""
    gen = {}
    order = []
    with open(g, 'r') as f:
        for line in f:
            tmp = line.strip().split()
            order.append(tmp[0])
            gen[tmp[0]] = tmp[1:]
    return (gen, order)


def compare_hets(iid, ped, geno):
    """Compare the genotypes of the parents and determine the expected number of
    segregating sites in the progeny. Then, count how many heterozygous sites
    there are in the progeny."""
    if iid not in ped:
        return ('NA', 'NA')
    p1 = ped[iid][0]
    p2 = ped[iid][1]
    if p1 not in geno or p2 not in geno:
        return ('NA', 'NA')
    if p1 == '0' and p2 == '0':
        exp_seg = 'NA'
    else:
        exp_seg = 0
        for site in zip(geno[p1], geno[p2]):
            # Skip sites with missing data
            if '9' in site:
                continue
            # If there's a het, then it should be segregating in the progeny
            elif '1' in site:
                exp_seg += 1
            # Last, if there is a homozygous difference between the parents,
            # then it is expected to be segregating
            elif site[0] != site[1]:
                exp_seg += 1
            else:
                continue
    # Then, count the observed heterozygous sites
    obs_het = geno[iid].count('1')
    return (exp_seg, obs_het)


def main(bopa, ped):
    """Main function."""
    pedigree = parse_ped(ped)
    genotypes, order = parse_geno(bopa)
    print 'ID NSeg NHet'
    for iid in order:
        e, o = compare_hets(iid, pedigree, genotypes)
        print iid, e, o
    return


if len(sys.argv) != 3:
    print """Summarize the pre-imputation genotype file that has the BOPA SNPs. We do it
this way because we want to keep track of the number of segregating sites that
are expected within each family. Takes two arguments:
    1) AlphaPeel BOPA input file
    2) Pedigree"""
    exit(1)
else:
    main(sys.argv[1], sys.argv[2])
