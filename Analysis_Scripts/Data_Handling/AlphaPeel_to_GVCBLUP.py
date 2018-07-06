#!/usr/bin/env python
"""Convert from the PED/MAP output from AlphaPeel to the input format needed
for GVCBLUP. Takes three arguments:
    1) Phenotyped lines file
    2) Single chromosome AlhapPeel PED (gzipped)
    3) Single chromosome AlhapPeel MAP (gzipped)
"""

import sys
import gzip


def t_ped(p, w):
    """Transpose the PED and return a list of lists with the genotypes"""
    n_ped = []
    samp = []
    with gzip.open(p, 'rt') as f:
        for line in f:
            tmp = line.strip().split()
            sname = tmp[1]
            genos = tmp[6:]
            if sname not in w:
                continue
            else:
                samp.append(sname)
                n_ped.append([a1+a2 for a1,a2 in zip(genos[::2], genos[1::2])])
    t_ped = list(zip(*n_ped))
    return (samp, t_ped)


def minor_allele(genos):
    """Figure out which allele is the minor allele at each site. Return them
    as a list."""
    mas = []
    maj = []
    for g in genos:
        calls = list(''.join(g))
        alleles = set(calls)
        alleles.discard('0')
        # Count the alleles
        a = list(alleles)
        if len(a) == 1:
            mas.append(a[0])
            maj.append('N')
        else:
            a1 = a[0]
            a2 = a[1]
            a1_c = calls.count(a1)
            a2_c = calls.count(a2)
            if a1_c >= a2_c:
                mas.append(a1)
                maj.append(a2)
            else:
                mas.append(a2)
                maj.append(a1)
    return (mas, maj)


def num_conv(mi, ma, tp):
    """Convert the transposed PED to 0/1/2 format using the minor alleles."""
    num = []
    for minor, major, genos in zip(mi, ma, tp):
        snp = []
        for ind in genos:
            if ind == minor + minor:
                snp.append('0')
            elif ind == minor + major or ind == major + minor:
                snp.append('1')
            elif ind == major + major:
                snp.append('2')
            else:
                snp.append('-9')
        num.append(snp)
    return num


def main(phen, ap_ped, ap_map):
    """Main function."""
    # Store the names of the phenotyped lines in a list
    with_phen = []
    with open(phen, 'r') as f:
        for line in f:
            with_phen.append(line.strip())
    with_phen = set(with_phen)
    # Parse the MAP file and store the SNP IDs in a list
    snp_ids = []
    with gzip.open(ap_map, 'rt') as f:
        for line in f:
            snp_ids.append(line.strip().split()[1])
    # Then parse the PED file.
    samples, tped = t_ped(ap_ped, with_phen)
    # Determine the minor alleles
    minor, major = minor_allele(tped)
    # Convert to 0/1/2
    numeric = num_conv(minor, major, tped)
    t_num = list(zip(*numeric))
    # Then print it back out
    print('Sample_ID ' + ' '.join(snp_ids))
    for s, g in zip(samples, t_num):
        print(' '.join([s] + list(g)))
    return


if len(sys.argv) != 4:
    print('''Convert from the PED/MAP output from AlphaPeel to the input format needed
for GVCBLUP. Takes three arguments:
    1) Phenotyped lines file
    2) Single chromosome AlhapPeel PED (gzipped)
    3) Single chromosome AlhapPeel MAP (gzipped)''')
    exit(1)
else:
    main(sys.argv[1], sys.argv[2], sys.argv[3])
