#!/usr/bin/env python
"""Script to trim down VCF based on intervals in a BED - bedtools SHOULD be able
to do this, but it fails on long references. Bedops similarly should be able to
do this, but it fails on VCFs with long lines (many samples). Weird. Writes
trimmed VCF data to stdout. Takes two arguments:
    1) VCF to be trimmed
    2) BED regions
"""

import sys
import pprint


def parse_bed(bedfile):
    """Parse the BED file and store it in a dictionary. The keys are the
    chromosome names and the values are lists of tuples for the regions."""
    bed_data = {}
    with open(bedfile, 'r') as f:
        for line in f:
            tmp = line.strip().split()
            interval = (int(tmp[1]), int(tmp[2]))
            if tmp[0] not in bed_data:
                bed_data[tmp[0]] = [interval]
            else:
                bed_data[tmp[0]].append(interval)
    return bed_data


def main(vcf, bed):
    """Main function"""
    bed_data = parse_bed(bed)
    #   Iterate through the VCF. Print header lines and print variants that
    #   fall in an interval.
    with open(vcf, 'r') as f:
        for line in f:
            if line.startswith('#'):
                print line.strip()
            else:
                tmp = line.strip().split('\t')
                chrom = tmp[0]
                pos = int(tmp[1])
                if chrom not in bed_data:
                    continue
                else:
                    for interval in bed_data[chrom]:
                        if pos > interval[0] and pos <= interval[1]:
                            print line.strip()
                            break


main(sys.argv[1], sys.argv[2])
