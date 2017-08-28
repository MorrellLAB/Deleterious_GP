#!/usr/bin/env python
"""Convert from Karen Beaubien's genotype calls to a VCF so that we can compare
the ALCHEMY calls to the curated calls. Takes three arguments:
    1) Karen's genotype output file (gzipped)
    2) BOPA physical VCF
    3) AB->ATCG file from Emily
"""

import sys
import gzip
import pprint


def fix_sample_name(name):
    """Fix the sample names from the KB report so that they match the sample
    names used in other files."""
    if name.startswith('G'):
        f1 = name.replace('G10W', 'MS10S3')
        fixed = f1.replace('-', '-0')
    elif name.startswith('MS11'):
        fixed = name
    elif name.startswith('MS12'):
        fixed = name
    return fixed

def parse_kb_calls(kb_call):
    """Parse Karen's genotype calls and store them as a dictionary:
    {
        Sample: {SNP_1: AA/BB/AB/00, SNP_2: AA/BB/AB/00, ...},
        ...
    }"""
    calls = {}
    with gzip.open(kb_call, 'rb') as f:
        for index, line in enumerate(f):
            if index < 9:
                continue
            elif index == 9:
                samples = line.strip().split()
                # fix the samples so they are consistent
                samples = [fix_sample_name(s) for s in samples]
                # Initialize the sub-dictionaries
                for s_id in samples:
                    calls[s_id] = {}
            else:
                tmp = line.strip().split()
                snp = tmp[0]
                # Then, build the dictionary
                for s_id, genotype in zip(samples, tmp[1:]):
                    calls[s_id].update({snp: genotype})
    return calls


def parse_vcf(bopa):
    """Parse the BOPA VCF and store the physical information for each SNP."""
    pass


def parse_ab(ab_trans):
    """Parse the AB->ATCG table and store them as a dictionary:
    {
        SNP: {'A': ATCG, 'B'; ATCG},
        ...
    }"""
    pass


def main(kb_calls, bopa, ab_trans):
    """Main function."""
    curated = parse_kb_calls(kb_calls)
    pprint.pprint(curated)
    return


if len(sys.argv) != 4:
    print """Convert from Karen Beaubien's genotype calls to a VCF so that we can compare
the ALCHEMY calls to the curated calls. Takes three arguments:
    1) Karen's genotype output file (gzipped)
    2) BOPA physical VCF
    3) AB->ATCG file from Emily"""
    exit(1)
else:
    main(sys.argv[1], sys.argv[2], sys.argv[3])
