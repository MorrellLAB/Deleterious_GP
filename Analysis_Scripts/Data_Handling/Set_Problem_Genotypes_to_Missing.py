#!/usr/bin/env python
"""Set site with genotypes that are not in an allowable list to missing calls.
Takes one argument:
    1) VCF to clean"""

import sys

#   The allowable genotypes
allowable = ['0/0', '0/1', '1/0', '1/1', './.']

with open(sys.argv[1], 'r') as f:
    for line in f:
        if line.startswith('#'):
            print line.strip()
        else:
            tmp = line.strip().split('\t')
            #   Iterate through the genotypes, if any genotype is not in the
            #   list of allowable genotypes, then we print the SNP ID and
            #   continue
            toprint = tmp[0:9]
            for g in tmp[9:]:
                if g == '.':
                    gt = './.'
                else:
                    geno = g.split(':')[0]
                    ad = g.split(':')[1]
                    if len(ad.split(',')) > 2:
                        ad_mod = '.'
                    else:
                        ad_mod = ad
                    if geno not in allowable:
                        gt = './.:' + ad_mod + ':' + ':'.join(g.split(':')[2:-1]) + ':.'
                    else:
                        t = g.split(':')
                        gt = ':'.join([t[0], ad_mod] + t[2:-1]) + ':.'
                toprint.append(gt)
            print '\t'.join(toprint)
