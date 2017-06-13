#!/usr/bin/env python
"""Count the genotypes in a PED file, and print out a summary of them. We will
check these against the BOPA alleles from ref-seq mapping. Takes three
arguments:
    1) BOPA alleles VCF
    2) PED to summarize
    3) MAP of marker order
"""

import sys


def marker_order(map_file):
    """Return a list giving ordered markers in the MAP file. Basically just
    slices out the second column."""
    ordered = []
    with open(map_file, 'r') as f:
        for line in f:
            ordered.append(line.strip().split()[1])
    return ordered


def bopa_alleles(vcf):
    """Return a dictionary of the BOPA SNPs with their alleles.
    {
        SNP_1: (ref, alt),
        SNP_2: (ref, alt),
        ...
    }
    """
    bopa = {}
    with open(vcf, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                tmp = line.strip().split()
                bopa[tmp[2]] = (tmp[3], tmp[4])
    return bopa


def transpose_ped(ped):
    """Return a PED file, transposed, and with the metadata columns removed.
    This is so we can iterate over markers isntead of individuals."""
    tped = []
    with open(ped, 'r') as f:
        for line in f:
            tmp = line.strip().split()
            tped.append(tmp[6:])
    return zip(*tped)


def count_genotypes(ped_column):
    """Count up the diploid genotypes, and return them in a tuple. They will
    look like this:
    (AA, TT, CC, GG, AT, AC, AG, TC, TG, GC, Missing).
    They are un-phased."""
    counts = {
        'AA': 0,
        'TT': 0,
        'CC': 0,
        'GG': 0,
        'AT': 0,
        'AC': 0,
        'AG': 0,
        'CT': 0,
        'GT': 0,
        'CG': 0,
        'Miss': 0
        }
    # re-transpose the column, to count genotypes. It's ugly, but we need to do
    # this to 
    for a1, a2 in zip(*ped_column):
        if a1 == '0' and a2 == '0':
            counts['Miss'] += 1
        else:
            call = ''.join(sorted([a1, a2]))
            counts[call] += 1
    geno_vec = (
        str(counts['AA']),
        str(counts['TT']),
        str(counts['CC']),
        str(counts['GG']),
        str(counts['AT']),
        str(counts['AC']),
        str(counts['AG']),
        str(counts['CT']),
        str(counts['GT']),
        str(counts['CG']),
        str(counts['Miss'])
    )
    return geno_vec


def check_problem(ref, alt, geno_vec):
    """Return a true/false for whether the genotypes are problematic. That would
    be if they segregate for alleles that are not present forward or reverse
    orientation of the physical map alleles."""
    # Make a string from the ref and alt for lookup to find the columns that are
    # supposed to be zero. Recall that the order of the genotype vector is
    # (AA, TT, CC, GG, AT, AC, AG, TC, TG, GC, Missing)
    exp = ''.join(sorted([ref, alt]))
    AA = 0
    TT = 1
    CC = 2
    GG = 3
    AT = 4
    AC = 5
    AG = 6
    TC = 7
    TG = 8
    GC = 9
    expected_zero_counts = {
        'AT': [CC, GG, AC, AG, TC, TG, GC],
        'AC': [AT, AG, TC, GC],
        'AG': [AT, AC, TG, GC],
        'CT': [AT, AC, TG, GC],
        'CG': [AA, TT, AC, AG, AT, TC, TG],
        'GT': [AT, AG, TC, GC]
        }
    # Then, get the counts for all the genotypes, and ask if they are 0
    counts = []
    for col in expected_zero_counts[exp]:
        if geno_vec[col] == '0':
            counts.append(True)
        else:
            counts.append(False)
    # If they are all True, we have no problem.
    if all(counts):
        return 'OK'
    else:
        return 'Error'


def main(vcf, plinkped, plinkmap):
    """Main function."""
    ordered_snps = marker_order(plinkmap)
    bopa = bopa_alleles(vcf)
    tped = transpose_ped(plinkped)
    # Print a header
    print 'SNP\tVCF_Ref\tVCF_Alt\tAA\tTT\tCC\tGG\tAT\tAC\tAG\tTC\tTG\tGC\tMissing\tCheck'
    for i, snp in enumerate(ordered_snps):
        geno_vec = count_genotypes(tped[2*i:2*(i+1)])
        check = check_problem(bopa[snp][0], bopa[snp][1], geno_vec)
        toprint = [snp, bopa[snp][0], bopa[snp][1]] + list(geno_vec) + [check]
        print '\t'.join(toprint)


main(sys.argv[1], sys.argv[2], sys.argv[3])
