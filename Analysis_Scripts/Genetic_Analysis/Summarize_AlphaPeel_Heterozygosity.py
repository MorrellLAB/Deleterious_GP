#!/usr/bin/env python
"""Summarize the .haps file from an AlphaPeel run in terms of heterozygosity.
We will also want to keep track of the number of segregating sites in any
parental combination. We will keep track of this with 0, 0.5, and 1 for AA, Aa,
and aa, respectively. Takes tws arguments:
    1) .haps output file
    2) Pedigree file
"""

import sys


def parse_pedigree(p):
    """Parse the pedigree. For each individual, store its parents, if it has
    any in the experiment."""
    ped = {}
    with open(p, 'r') as f:
        for line in f:
            tmp = line.strip().split()
            ped[tmp[0]] = (tmp[1], tmp[2])
    return ped


def parse_haps(h):
    """Parse the haps and store the most likely state for each site in each
    individual. Because the haps file has four lines per individual, we have
    to iterate through the file four lines at a time."""
    haps = {}
    order = []
    with open(h, 'r') as f:
        for p_AA, p_Aa, p_aA, p_aa in zip(f, f, f, f):
            for index, site_probs in enumerate(zip(p_AA.split(), p_Aa.split(), p_aA.split(), p_aa.split())):
                # The first element is the sample ID.
                if index == 0:
                    sample = site_probs[0]
                    haps[sample] = []
                    order.append(sample)
                else:
                    hom_ref = float(site_probs[0])
                    het = float(site_probs[1]) + float(site_probs[2])
                    hom_alt = float(site_probs[3])
                    vals = [0.0, 0.5, 1.0]
                    probs = [hom_ref, het, hom_alt]
                    highest = max(probs)
                    haps[sample].append(vals[probs.index(highest)])
    return (haps, order)


def compare_hets(iid, pedi, genotypes):
    """Compare the genotypes of the progeny to the genotypes of the parents and
    determine the number of segregating sites in the parents, and how many
    heterozygous sites there are in the progeny. Because the progeny are F3s, we
    would expect about 12.5% heterozygosity."""
    # get the parents
    if iid not in pedi:
        return ('NA', 'NA')
    par1 = pedi[iid][0]
    par2 = pedi[iid][1]
    # If the individual is a founder (both parents 0), then return NAs for the
    # number of seg sites
    if par1 == '0' and par2 == '0':
        exp_seg = 'NA'
    else:
        exp_seg = 0
        for site in zip(genotypes[par1], genotypes[par2]):
            # We will take advantage of the integer type here:
            # if p1+p2 != 0 and 2, then it will be segregating
            x = site[0] + site[1]
            if x != 0 and x != 2:
                exp_seg += 1
    # Then, count the actual number of heterozygous sites. This is just the
    # count of 0.5 values in the individual
    obs_het = genotypes[iid].count(0.5)
    return (exp_seg, obs_het)


def main(haps, ped):
    """Main function."""
    pedigree = parse_pedigree(ped)
    # Then parse the haps file
    genotypes, order = parse_haps(haps)
    # Then, for each individual, in order, print the expected number of
    # segregating sites and the o bserved number of heterzoygous sites
    print 'ID NSeg NHet'
    for iid in order:
        e, o = compare_hets(iid, pedigree, genotypes)
        print iid, e, o
    return


if len(sys.argv) != 3:
    print """Summarize the .haps file from an AlphaPeel run in terms of heterozygosity.
We will also want to keep track of the number of segregating sites in any
parental combination. We will keep track of this with 0, 0.5, and 1 for AA, Aa,
and aa, respectively. Takes tws arguments:
    1) .haps output file
    2) Pedigree file"""
    exit(1)
else:
    main(sys.argv[1], sys.argv[2])
