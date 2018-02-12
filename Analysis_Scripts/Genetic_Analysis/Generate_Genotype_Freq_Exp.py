#!/usr/bin/env python
"""Generate an expectation for genotype frequencies, given the pedigree and
parental genotypes for each marker. We will also take into account that there
are two rounds of inbreeding following the initial crosses. Takes four
arguments:
    1) PED file with genotypes
    2) Pedigree file from AlphaPeel
    3) PLINK BIM file for allele assignment
    4) Cycle (1 for C1 expectaion, 2 for C2 expectation, 3 for C3)
"""

import sys


def parse_pedigree(p, c):
    """Parse the pedigree and return a list of the indivudal IDs that are used
    in the crosses."""
    # Check the cycle that was passed.
    if c == '1':
        k = 'MS10'
    elif c == '2':
        k = 'MS11'
    elif c == '3':
        k = 'MS12'
    else:
        sys.stderr.write('Please supply 1, 2, or 3 for the cycle.\n')
        exit(2)
    matings = []
    with open(p, 'r') as f:
        for line in f:
            if not line.startswith(k):
                continue
            else:
                tmp = line.strip().split()
                p1 = tmp[1]
                p2 = tmp[2]
                if (p1, p2) in matings:
                    continue
                else:
                    matings.append((p1, p2))
    return matings


def extract_geno(ped, ids):
    """Return a dictionary of genotypes for each parent given."""
    geno = {}
    with open(ped, 'r') as f:
        for line in f:
            tmp = line.strip().split()
            if tmp[1] not in ids:
                continue
            else:
                snps = tmp[6:]
                nsnps = len(snps)/2
                g = [(snps[2*i], snps[(2*i)+1]) for i in range(nsnps)]
                geno[tmp[1]] = g
    return geno


def read_bim(b):
    """Read the bim file from the PLINK binary output. A1 is the 5th field, and
    A2 is the 6th field."""
    a = []
    names = []
    with open(b, 'r') as f:
        for line in f:
            tmp = line.strip().split()
            a.append((tmp[4], tmp[5]))
            names.append(tmp[1])
    return a, names


def expected(c, p, a, famsize=24):
    """For each cross, generate the number of expected genotypes in each class,
    then return a sum across all families for the cycle of interest. For this,
    we are assuming an F3 family, because it is difficult to generalize the
    expected genotype numbers on-the-fly."""
    by_fam = []
    for cross in c:
        f = []
        p1 = p[cross[0]]
        p2 = p[cross[1]]
        for p1_g, p2_g, alleles in zip(p1, p2, a):
            # for ease, convert the parental calls into A/B
            if p1_g == (alleles[0], alleles[0]):
                p1_ab = 'AA'
            elif p1_g == (alleles[0], alleles[1]) or p1_g == (alleles[1], alleles[0]):
                p1_ab = 'AB'
            elif p1_g == (alleles[1], alleles[1]):
                p1_ab = 'BB'
            # The same with P2
            if p2_g == (alleles[0], alleles[0]):
                p2_ab = 'AA'
            elif p2_g == (alleles[0], alleles[1]) or p2_g == (alleles[1], alleles[0]):
                p2_ab = 'AB'
            elif p2_g == (alleles[1], alleles[1]):
                p2_ab = 'BB'
            # Then, enumerate over the possible genotpyes. This is annoying.
            if p1_ab == 'AA' and p2_ab == 'AA':
                aa = 1.0 * famsize
                ab = 0
                bb = 0
            elif (p1_ab == 'AA' and p2_ab == 'AB') or ( p1_ab == 'AB' and p2_ab == 'AA'):
                aa = (11.0/16.0) * famsize
                ab = (2.0/16.0) * famsize
                bb = (3.0/16.0) * famsize
            elif p1_ab == 'AB' and p2_ab == 'AB':
                aa = (7.0/16.0) * famsize
                ab = (2.0/16.0) * famsize
                bb = (7.0/16.0) * famsize
            elif (p1_ab == 'AA' and p2_ab == 'BB') or (p1_ab == 'BB' and p2_ab == 'AA'):
                aa = (3.0/8.0) * famsize
                ab = (2.0/8.0) * famsize
                bb = (3.0/8.0) * famsize
            elif (p1_ab == 'BB' and p2_ab == 'AB') or ( p1_ab == 'AB' and p2_ab == 'BB'):
                bb = (11.0/16.0) * famsize
                ab = (2.0/16.0) * famsize
                aa = (3.0/16.0) * famsize
            elif p1_ab == 'BB' and p2_ab == 'BB':
                bb = 1.0 * famsize
                ab = 0
                aa = 0
            # Then put them into the currnet family list
            f.append((aa, ab, bb))
        # At this point, f is a list of lists, with the sublists representing
        # the expected numbers of aa, ab, and bb progeny
        by_fam.append(f)
    # Then, sum across all families
    cycle_tot = []
    for fams in zip(*by_fam):
        aa_tot = sum([x[0] for x in fams])
        ab_tot = sum([x[1] for x in fams])
        bb_tot = sum([x[2] for x in fams])
        cycle_tot.append((aa_tot, ab_tot, bb_tot))
    return cycle_tot


def main(geno, pedigree, bim, cycle):
    """Main function."""
    # First, parse the pedigree so that we know which sections of the genotype
    # matrix we need to read. We will only keep the diploid genotypes of the
    # parental lines.
    crosses = parse_pedigree(pedigree, cycle)
    # Then, "flatten" the list of crosses to get the parental IDs
    parents = list(set([par for cross in crosses for par in cross]))
    # Next, parse the genotypes for these parents
    p_geno = extract_geno(geno, parents)
    # Read the alleles from the bim file
    alleles, snp_names = read_bim(bim)
    # Then, generate the expected genotype counts for each family
    exp_geno_counts = expected(crosses, p_geno, alleles)
    # Fianlly! Print out the expected counts
    print 'SNPName ExpAA ExpAB ExpBB'
    for n, c in zip(snp_names, exp_geno_counts):
        print n, c[0], c[1], c[2]
    return


if len(sys.argv) != 5:
    print """Generate an expectation for genotype frequencies, given the pedigree and
parental genotypes for each marker. We will also take into account that there
are two rounds of inbreeding following the initial crosses. Takes four
arguments:
    1) PED file with genotypes
    2) Pedigree file from AlphaPeel
    3) PLINK BIM file for allele assignment
    4) Cycle (1 for C1 expectaion, 2 for C2 expectation, 3 for C3)"""
    exit(1)
else:
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
