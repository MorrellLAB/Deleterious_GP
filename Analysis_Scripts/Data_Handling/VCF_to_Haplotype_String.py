#!/usr/bin/env python
"""Convert from a phased VCF to a ATCG phased haplotype string. We are doing
this because we want to be able to easily check the haplotype output from
beagle against that of PHASE. Takes one argument:
    1) Phased VCF
"""

import sys


def main(vcf):
    """Main function."""
    with open(vcf, 'r') as f:
        for line in f:
            if line.startswith('##'):
                continue
            elif line.startswith('#CHROM'):
                samples = line.strip().split()[9:]
                haplotypes = {}
                for s in samples:
                    haplotypes[s] = [[], []]
            else:
                tmp = line.strip().split()
                ref = tmp[3]
                alt = tmp[4]
                # We assume that the GT field is the first in the sample-level
                # fields
                gts = [x.split(':')[0] for x in tmp[9:]]
                # Start counting up the haploypes
                for s, g in zip(samples, gts):
                    if '|' not in g:
                        sys.stderr.write('This VCF does not appear to be phased!\n')
                        exit(1)
                    else:
                        h = g.split('|')
                        for i, allele in enumerate(h):
                            if allele == '0':
                                haplotypes[s][i].append(ref)
                            elif allele == '1':
                                haplotypes[s][i].append(alt)
                            else:
                                haplotypes[s][i].append('?')
    # Then, print it out
    for samp in sorted(haplotypes):
        print samp
        print ' '.join(haplotypes[samp][0])
        print ' '.join(haplotypes[samp][1])
    return


if len(sys.argv) != 2:
    print """Convert from a phased VCF to a ATCG phased haplotype string. We are doing
this because we want to be able to easily check the haplotype output from
beagle against that of PHASE. Takes one argument:
    1) Phased VCF"""
    exit(1)
else:
    main(sys.argv[1])
