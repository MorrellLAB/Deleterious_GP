#!/usr/bin/env python
"""Convert PHASE output to a phased reference panel for BEAGLE imputation. This
script will use the "pairs" output from PHASE to do this, and assume that the
most likely pair is the "true" configuration. NOTE, this script assumes that the
PHASE SNPs and the BOPA SNPs in the physical VCF are in the same order. Takes
three arguments:
    1) PHASE .pairs file
    2) BOPA physical VCF
    3) Chromosome name
"""

import sys
import pprint


def parse_pairs(phasepairs):
    """Parse the PHASE pairs and return the ordered haplotypes for each
    individual. Will just take the most likely pairs (first listed)."""
    p = {}
    with open(phasepairs, 'r') as f:
        for line in f:
            if line.startswith('IND:'):
                # They are coded as IND: #ID, so we strip off the first six
                # characters
                ind_id = line.strip()[6:]
                # Then, read the next line, which will have the haplotypes
                haps = f.next().split()
                # The haplotypes are separated by a space-comma-space, so the
                # first and the third fields have the haplotypes. It is given
                # as one string for each chromosome, but we want to make it
                # such that the genotypes are together. First element will be
                # haplotype 1, and the second element will be haplotype 2.
                p[ind_id] = []
                for h1, h2 in zip(haps[0], haps[2]):
                    p[ind_id].append((h1, h2))
    return p


def parse_vcf(vcf, chrom):
    """Parse the VCF, and only save SNPs from the specified chromosome. We want
    to keep all the VCF metadata."""
    positions = []
    with open(vcf, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                tmp = line.strip().split()
                # Check the chromosome. If we are on the wrong chromosome, then
                # just skip the SNP
                c = tmp[0]
                if chrom == c:
                    # Save the physical position, ref, alt, and ID
                    metadata = (tmp[1], tmp[2], tmp[3], tmp[4])
                    positions.append(metadata)
                else:
                    continue
    return positions


def generate_vcf(haplotypes, chrom, vcfinfo):
    """Use the metadata from the physical VCF and the phased haplotypes to
    generate a phased VCF to use as a reference panel. Print it to stdout."""
    samples = sorted(haplotypes)
    # Define the VCF header
    header = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t""" + '\t'.join(samples)
    print header
    for index, snp in enumerate(vcfinfo):
        # keep the meta data 
        toprint = [chrom, snp[0], snp[1], snp[2], snp[3], '.', '.', '.', 'GT']
        # Then start printing the haplotypes
        for sample in sorted(haplotypes):
            gt = haplotypes[sample][index]
            h1 = gt[0]
            h2 = gt[1]
            if h1 == snp[2]:
                h1 = '0'
            elif h1 == snp[3]:
                h1 = '1'
            if h2 == snp[2]:
                h2 = '0'
            elif h2 == snp[3]:
                h2 = '1'
            toprint.append(h1 + '|' + h2)
        print '\t'.join(toprint)
    return


def main(pairs, bopa, chrom):
    """Main function."""
    # Parse the pairs information
    haps = parse_pairs(pairs)
    # Then parse the VCF
    physical = parse_vcf(bopa, chrom)
    # Generate and print the VCF
    generate_vcf(haps, chrom, physical)
    return


if len(sys.argv) != 4:
    print """Convert PHASE output to a phased reference panel for BEAGLE imputation. This
script will use the "pairs" output from PHASE to do this, and assume that the
most likely pair is the "true" configuration. NOTE, this script assumes that the
PHASE SNPs and the BOPA SNPs in the physical VCF are in the same order. Takes
three arguments:
    1) PHASE .pairs file
    2) BOPA physical VCF
    3) Chromosome name"""
    exit(1)
else:
    main(sys.argv[1], sys.argv[2], sys.argv[3])
