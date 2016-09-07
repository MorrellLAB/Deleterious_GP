#!/usr/bin/env python
"""Extract the nonsynonymous variants from ANNOVAR output, and create
substittions files for PolyPhen2, PROVEAN, and BAD_Mutations."""

import sys

ns_snps = {}
with open(sys.argv[1], 'r') as f:
    for line in f:
        tmp = line.strip().split('\t')
        impact = tmp[1]
        if 'nonsynonymous' in impact:
            ann = tmp[2].split(',')[0].split(':')
            txid = ann[1]
            subs = ann[4]
            aa1 = subs[2]
            aa2 = subs[-1]
            pos = subs[3:-1]
            if txid not in ns_snps:
                ns_snps[txid] = [(aa1, aa2, pos)]
            else:
                ns_snps[txid].append((aa1, aa2, pos))

#   Make the PPH2 subs
handle = open('PPH2_Subs.txt', 'w')
for k, v in ns_snps.iteritems():
    snps = sorted(v, key=lambda x: int(x[2]))
    for entry in snps:
        handle.write('\t'.join([k] + list(entry)) + '\n')
handle.close()

#   And make the PROVEAN subs
for k, v in ns_snps.iteritems():
    handle = open('PROVEAN/' + k + '.subs', 'w')
    snps = sorted(v, key=lambda x: int(x[2]))
    for entry in snps:
        handle.write(entry[0] + entry[2] + entry[1] + '\n')
    handle.close()

#   Make the BAD_Mutations subs
for k, v in ns_snps.iteritems():
    handle = open('DM/' + k + '.subs', 'w')
    snps = sorted(v, key=lambda x: int(x[2]))
    for entry in snps:
        handle.write(entry[2] + '\n')
    handle.close()
