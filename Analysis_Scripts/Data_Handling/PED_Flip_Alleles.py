#!/usr/bin/env python
"""Parse the PED file for the genomic prediction experiment, and flip the
alleles so that they accord with the forward strand of the reference, as
determined by mapping of BOPA contextual sequences to the Morex assembly. Takes
three arguments:
    1) PED file
    2) MAP file
    3) BOPA vcf
"""


import sys


def marker_order(mapfile):
    """Return a list of ordered markers."""
    ordered = []
    with open(mapfile, 'r') as f:
        for line in f:
            ordered.append(line.strip().split()[1])
    return ordered


def bopa_alleles(vcf):
    """Return a dictionary of BOPA SNPs and their alleles."""
    bopa = {}
    with open(vcf, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                tmp = line.strip().split()
                bopa[tmp[2]] = (tmp[3], tmp[4])
    return bopa


def parse_ped(ped, parents=21):
    """Return pieces of a PED file for parsing. Will return the metadata,
    the parents of the population, and the progeny in separate pieces."""
    meta = []
    t_par = []
    t_prog = []
    with open(ped, 'r') as f:
        for index, line in enumerate(f):
            tmp = line.strip().split()
            meta.append(tmp[:6])
            if index < parents:
                t_par.append(tmp[6:])
            else:
                t_prog.append(tmp[6:])
    # Return the transposed genotype matrices so that we can iterate over
    # markers rather than individuals
    return meta, zip(*t_par), zip(*t_prog)


def check_alleles(alleles, observed):
    """Check the observed genotypes against the alleles from the VCF. We should
    only observe four genotype states: hom ref, hom alt, het, missing. If there
    are others, then we will reverse-complement the genotypes."""
    revcomp = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C'
        }
    ex = set([
        alleles[0] + alleles[0],
        alleles[0] + alleles[1],
        alleles[1] + alleles[0],
        alleles[1] + alleles[1],
        '00'])
    # Get the list of observed genotypes
    ob = set([
        g[0] + g[1]
        for g
        in zip(*observed)
        ])
    print alleles, ex, ob, ob.issubset(ex)
    return


def main(pedfile, mapfile, vcf):
    """Main function."""
    ordered_snps = marker_order(mapfile)
    bopa = bopa_alleles(vcf)
    meta, parents, progeny = parse_ped(pedfile)
    # For each genotype column, check it against the BOPA alleles
    for i, snp in enumerate(ordered_snps):
        check_alleles(bopa[snp], parents[2*i:2*(i+1)])
        check_alleles(bopa[snp], progeny[2*i:2*(i+1)])
    return

# def flipstrand(alleles, pedcols):
#     """Reverse-complement the alleles in the PED column"""
#     # Define a reverse complementation dictionary for flipping
#     revcomp = {
#         'A': 'T',
#         'T': 'A',
#         'G': 'C',
#         'C': 'G'
#         }
#     fixed_alleles = []
#     # Iterate through the genotype calls, and flip them
#     for geno in zip(*pedcols):
#         if geno[0] == '0' and geno[1] == '0':
#             fixed_alleles.append(('0', '0'))
#         elif geno == (revcomp[alleles[0]], revcomp[alleles[0]]):
#             # If the observed genotype is the reverse complement of the ref
#             # forward strand allele, flip it to the forward ref allele.
#             fixed_alleles.append((alleles[0], alleles[0]))
#         elif geno == (revcomp[alleles[1]], revcomp[alleles[1]]):
#             # If the observed genotype is the reverse complement of the forward
#             # alternate allele, set it to the forward alternate allele
#             fixed_alleles.append((alleles[1], alleles[1]))
#         elif geno == (revcomp[alleles[0]], revcomp[alleles[1]]) or geno == (revcomp[alleles[1]], revcomp[alleles[0]]):
#             # If the genotype is a heterozygous reverse complement, set it
#             # to the heterozygous forward state
#             fixed_alleles.append((alleles[0], alleles[1]))
#         else:
#             # Otherwise, just leave it
#             fixed_alleles.append(geno)
#     return zip(*fixed_alleles)

main(sys.argv[1], sys.argv[2], sys.argv[3])
