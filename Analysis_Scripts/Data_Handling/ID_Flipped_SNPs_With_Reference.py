#!/usr/bin/env python
"""Script to identify SNPs from Morex genotyping data that are RCed with
respect to the reference genome. Assumes that the Morex calls are homozygous
only. Takes two arguments:
    1) Morex genotyping data
    2) VCF of the SNPs against the reference
"""

import sys


def parse_geno(genotypes):
    """Store the genotyping data in a dictionary of the form
    {
        SNP_ID: haploid genotype,
        ...
    }
    """
    geno = {}
    with open(genotypes, 'r') as f:
        for index, line in enumerate(f):
            if index == 0:
                continue
            else:
                tmp = line.strip().split()
                geno[tmp[0]] = tmp[1]
    return geno


def parse_vcf(vcf):
    """Store the VCF reference allele in a dictionary of the form
    {
        SNP_ID: ref,
        ...
    }
    """
    vcf_alleles = {}
    with open(vcf, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                tmp = line.strip().split()
                vcf_alleles[tmp[2]] = tmp[3]
    return vcf_alleles


def compare_calls(genotype, vcf):
    """Compare the genotype call to the reference genome call. Returns one of
    three values:
        Match: Morex genotype matches the reference genome forward strand
        Flip: Morex genotype is the RC of the reference genome
        Problem: Morex genotype is not forward, nor RC
    """
    complement = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C'
        }
    if genotype == vcf:
        return "Match"
    elif genotype == complement[vcf]:
        return "Flip"
    else:
        return "Problem"


def main(genotyping, vcf):
    """Main function."""
    genotypes = parse_geno(genotyping)
    vcf_states = parse_vcf(vcf)
    # Print a header
    print 'SNP_ID\tGenotype_Call\tReference_Call\tStatus'
    for snp in sorted(vcf_states):
        if snp in genotypes:
            # Check for missing
            if genotypes[snp] == 'N':
                s = 'NA'
            else:
                s = compare_calls(genotypes[snp], vcf_states[snp])
            print '\t'.join([snp, genotypes[snp], vcf_states[snp], s])
    return


if len(sys.argv) != 3:
    print """
    Script to identify SNPs from Morex genotyping data that are RCed with
respect to the reference genome. Assumes that the Morex calls are homozygous
only. Takes two arguments:
    1) Morex genotyping data
    2) VCF of the SNPs against the reference"""
    exit(1)
else:
    main(sys.argv[1], sys.argv[2])
