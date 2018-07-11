#!/usr/bin/env python
"""Trim down the AlphaPeel dosages file by a given list of SNPs. Takes four
arguments:
    1) List of SNPs to keep (Gzipped)
    2) List of lines to keep (PLINK format)
    3) AlphaPeel SNP list
    4) (Gzipped) AlphaPeel dosage file
"""

import sys
import gzip


def main(snp_keep, line_keep, ap_snps, ap_dosages):
    """Main function."""
    # Read in a list of SNPs to keep, save it as a set
    k = []
    with gzip.open(snp_keep, 'rt') as f:
        for line in f:
            k.append(line.strip())
    ks = set(k)
    l = []
    with open(line_keep, 'r') as f:
        for line in f:
            l.append(line.strip().split()[1])
    l = set(l)
    # Make a list of offsets to keep
    k_i = []
    k_n = []
    with open(ap_snps, 'r') as f:
        for index, line in enumerate(f):
            tmp = line.strip().split()
            if tmp[2] in ks:
                k_i.append(index)
                k_n.append(tmp[2])
    k_i = set(k_i)
    # Then, start iterating through the file
    sys.stdout.write(' '.join(['IID'] + k_n) + '\n')
    with gzip.open(ap_dosages, 'rt') as f:
        for line in f:
            tmp = line.strip().split()
            if tmp[0] not in l:
                continue
            else:
                t = [tmp[0]]
                for index, v in enumerate(tmp[1:]):
                    if index in k_i:
                        t.append(v)
                sys.stdout.write(' '.join(t) + '\n')
                sys.stdout.flush()
    return



if len(sys.argv) != 5:
    print("""Trim down the AlphaPeel dosages file by a given list of SNPs. Takes three
arguments:
    1) List of SNPs to keep
    2) List of lines to keep
    3) AlphaPeel SNP list
    4) (Gzipped) AlphaPeel dosage file""")
    exit(1)
else:
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
