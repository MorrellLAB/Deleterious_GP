#!/usr/bin/env python
"""Read a VCF and calculate the proportion of homozygous major, heterozygous,
and homozygous minor genotypes for each minor allele count class."""

import sys

mac_genos = {}
with open(sys.argv[1], 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        else:
            tmp = line.strip().split('\t')
            if len(tmp[3]) != 1 or len(tmp[4]) != 1:
                continue
            else:
                gt = [g.split(':')[0].split('/') for g in tmp[9:]]
                alleles = [a for c in gt for a in c]
                #   Calculate the minor allele
                numref = alleles.count('0')
                numalt = alleles.count('1')
                if numref > numalt:
                    minor = '1'
                    major = '0'
                    mac = numalt
                elif numalt < numref:
                    minor = '0'
                    major = '1'
                    mac = numref
                elif numref == numalt:
                    minor = '1'
                    major = '0'
                    mac = numref
                #   Then ask if the minor allele count is present in the dict
                if mac not in mac_genos:
                    mac_genos[mac] = {
                        'HomMinor': 0,
                        'Het': 0,
                        'HomMajor': 0}
                for g in gt:
                    if g == [major, major]:
                        mac_genos[mac]['HomMajor'] += 1
                    if g == [minor, minor]:
                        mac_genos[mac]['HomMinor'] += 1
                    if g == [major, minor] or g == [minor, major]:
                        mac_genos[mac]['Het'] += 1

print 'MAC\tPHomMin\tPHet\tPHomMaj'
for i in sorted(mac_genos):
    totalgeno = float(mac_genos[i]['HomMajor'] + mac_genos[i]['Het'] + mac_genos[i]['HomMinor'])
    print str(i) + '\t' + str(mac_genos[i]['HomMinor']/totalgeno) + '\t' + str(mac_genos[i]['Het']/totalgeno) + '\t' + str(mac_genos[i]['HomMajor']/totalgeno)
